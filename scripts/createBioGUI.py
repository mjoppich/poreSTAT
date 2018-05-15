import argparse
import re

import os
import string

import sys
from enum import Enum

from lxml import etree

import yaml
from yaml import scanner
from jinja2 import Environment
from jinja2.loaders import FileSystemLoader

from porestat.utils.ArgParseExt import FolderType, FileStubType

import sys
from yaml import load
from collections import OrderedDict


class InputBinding:
    def __init__(self, ib):
        self.position = ib.get('position', None)
        self.prefix = ib.get('prefix', None)


class OutputBinding:
    def __init__(self, ob):
        self.glob = ob.get('glob', None)


class Param:
    optional = False
    default = None
    type = None

    def get_type(self):
        return self.type


class InputParam(Param):
    def __init__(self, param):
        self.id = param['id']
        self.type = param.get('type', None)
        if type(self.type) is str and self.type[-2:] == '[]': # v.1.0 syntax simplification ('<type>[]' == array of <type> )
            self.type = "array"
        if (type(self.type) is list and self.type[0] == 'null'):
            self.optional = True
        elif type(self.type) is str and self.type[-1] == '?':  # v.1.0 ('<type>?' == ['null', '<type>'])
            self.optional = True
            self.type = self.type[:-1]
        else:
            self.optional = False
        self.description = param.get('doc', param.get('description', None))
        self.default = param.get('default', None)
        input_binding = param.get('inputBinding', None)
        if input_binding:
            self.input_binding = InputBinding(input_binding)

    def get_type(self):
        if type(self.type) is list and self.type[0] == 'null':
            arg_type = self.type[1]
            if type(arg_type) is dict:
                return arg_type['type']
            else:
                return arg_type
        elif type(self.type) is dict:
            return self.type['type']
        else:
            return self.type

    def is_required(self):
        if type(self.type) is list and self.type[0] == 'null':
            return False
        return True

class OutputParam(Param):
    def __init__(self, param):
        self.id = param['id']
        self.type = param.get('type', None)
        self.description = param.get('description', None)
        output_binding = param.get('outputBinding', None)
        if output_binding:
            self.output_binding = OutputBinding(output_binding)


class Tool:
    def __init__(self, filename):

        with open(filename) as f:
            tool = load(f)
        try:
            self.tool_class = tool['class']
        except KeyError:
            sys.exit('`class` attribute of the CWL document not found')
        if self.tool_class != 'CommandLineTool':
            raise ValueError('Wrong tool class')
        try:
            self.basecommand = tool['baseCommand']
        except KeyError:
            sys.exit('`baseCommand` attribute of the CWL document not found')
        self.inputs = OrderedDict()
        if type(tool['inputs']) is list:  # ids not mapped
            for param_dict in tool['inputs']:
                param = InputParam(param_dict)
                self.inputs[param.id] = param
        elif type(tool['inputs']) is dict:  # ids mapped
            for id, param_dict in tool['inputs'].items():
                param_dict['id'] = id
                param = InputParam(param_dict)
                self.inputs[id] = param

        self.outputs = OrderedDict()
        if tool['outputs']:
            if type(tool['outputs']) is list:  # ids not mapped
                for param_dict in tool['outputs']:
                    param = OutputParam(param_dict)
                    self.outputs[param.id] = param
            elif type(tool['outputs']) is dict:  # ids mapped
                for id, param_dict in tool['outputs'].items():
                    param_dict['id'] = id
                    param = OutputParam(param_dict)
                    self.outputs[id] = param
        self.description = tool.get('doc', tool.get('description', None))
        self.cwl_version = tool.get('cwlVersion', '')


class Argument:
    def __init__(self, arg, fileOp=None):
        self.dest = Argument._get_dest(arg)
        self.help = arg.description
        self.option_strings = Argument._check_conflicting_prefixes(Argument._get_option_string(arg))
        self.default = Argument._get_default(arg)
        self.action = Argument._get_actions(arg)
        self.type = Argument._get_type(arg, fileOp)
        self.nargs = Argument._get_nargs(arg)
        self.choices = Argument._get_choices(arg)
        self.required = Argument._get_required(arg)
        if self.action:
            self.type = self.default = None

    @staticmethod
    def _get_dest(arg):
        s = arg.id.strip(string.punctuation)
        if arg.prefix:
            s = arg.prefix + s
        return re.sub(r'[{0}]'.format(string.punctuation), '_', s)

    @staticmethod
    def _get_option_string(arg):
        if arg.optional:
            if hasattr(arg, 'input_binding'):
                if arg.input_binding.prefix:
                    name = arg.input_binding.prefix.strip(string.punctuation)
                    if len(name) == 1:
                        return '-' + name
                    else:
                        return '--' + name
                else:
                    return '--' + arg.id.strip(string.punctuation)
        else:
            return Argument._get_dest(arg)

    @staticmethod
    def _check_conflicting_prefixes(name):
        global argument_names
        if name in argument_names:
            # if the name already exists, add '_N' to it, where N is an order number of this name
            # example: if argument 'foo' already exists, an argument 'foo_1' is created
            # if 'foo_1' exists, 'foo_2' is created
            same_names = list(filter(lambda x: x.startswith(name), argument_names))
            return name + '_{0}'.format(len(same_names))
        else:
            argument_names.append(name)
            return name


    @staticmethod
    def _get_required(arg):
        return arg.is_required()

    @staticmethod
    def _get_type(arg, fileOp):
        CWL_TO_PY_TYPES = {
            'string': 'str',
            'int': 'int',
            'boolean': 'bool',
            'double': 'float',
            'float': 'float',
            'array': 'list',
            'File': 'argparse.FileType({fileOp})'.format(fileOp=fileOp),
            'stdout': 'argparse.FileType({fileOp})'.format(fileOp=fileOp),
            'stderr': 'argparse.FileType({fileOp})'.format(fileOp=fileOp),
            'enum': 'str',
            'Directory': 'FolderType({fileOp})'.format(fileOp=fileOp)
        }

        argtype = None
        required=True

        if isinstance(arg.type, str):
            argtype = CWL_TO_PY_TYPES[arg.type]

        elif isinstance(arg.type, list):

            if arg.type.index('Null') == 0 and len(arg.type.index) > 1:
                argtype = CWL_TO_PY_TYPES[arg.type[1]]


        elif isinstance(arg.type, dict):
            itemtype = arg.type['items']
            argtype = CWL_TO_PY_TYPES[itemtype]

        if argtype == None:

            arg_type = CWL_TO_PY_TYPES[arg.get_type()]
            if arg_type is list and type(arg_type) is list:
                return None
            else:
                return arg_type

        else:
            return argtype

    @staticmethod
    def _get_choices(arg):
        if arg.get_type() == 'enum':
            if type(arg.type) is list and arg.type[0] == 'null':
                return arg.type[1]['symbols']
            elif type(arg.type) is dict:
                return arg.type['symbols']

    @staticmethod
    def _get_default(arg):
        if arg.default:
            if type(arg.default) is str:
                return "\"" + arg.default + "\""  # for proper rendering in j2 template
            else:
                return arg.default

    @staticmethod
    def _get_actions(arg):
        if arg.optional and arg.get_type() == 'boolean':
            return 'store_true'

    @staticmethod
    def _get_nargs(arg):
        if arg.type == 'array':
            if arg.optional:
                return '*'
            else:
                return '+'

class GUIObjects(Enum):
    FILEDIALOG=0
    LABEL=1
    BUTTON=2
    GROUPBOX=3
    TEXTINPUT=4
    CHECKBOX=5


class ActionElem:

    def __init__(self, elem, prefix="parser1"):

        self.original_elem = elem
        self.prefix = prefix


    def get_nargs(self):

        if self.original_elem.nargs != None:
            return self.original_elem.nargs

        if isinstance(self.original_elem.type, (str, Enum, int)):
            return 1

        #default value
        return 1

    def get_type(self):

        if isinstance(self.original_elem.type, (int)):
            return int
        elif isinstance(self.original_elem.type, (float)):
            return float
        elif isinstance(self.original_elem.type, (str)):
            return str

        return None

    def get_name(self):

        optStr = ""

        for x in self.original_elem.option_strings:
            if len(x) > len(optStr):
                optStr = x

        return optStr

    def get_gui_object_type(self):

        if isinstance(self.original_elem.type, (argparse.FileType, FileStubType, FolderType)):
            return GUIObjects.FILEDIALOG

        if self.original_elem.type == None and self.get_nargs() == 0:
            return GUIObjects.CHECKBOX

        return GUIObjects.TEXTINPUT

    def make_biogui_component(self):

        componentType = self.get_gui_object_type()

        # case 1: simple flag, e.g. --help
        if componentType == GUIObjects.FILEDIALOG:
            return self.make_filedialog()
        elif componentType == GUIObjects.CHECKBOX:
            return self.make_checkbox()
        elif componentType == GUIObjects.TEXTINPUT:
            return self.make_text_input()

        return (None, None, None)

    def get_help(self):
        return self.original_elem.help if self.original_elem.help != None else ""

    def get_dest(self):
        return self.original_elem.dest

    def get_default(self):
        return self.original_elem.default

    def is_required(self):
        return self.original_elem.required if self.original_elem.required != None else False

    def make_filedialog(self):

        groupElem = etree.Element('group')
        baseElem = etree.SubElement(groupElem, "hgroup")

        labelElem = etree.SubElement(baseElem, "label")
        labelElem.text = self.get_dest()

        outputMode = 'false'

        if 'w' in self.original_elem.type._mode:
            outputMode = 'true'

        folderMode = 'false'
        if isinstance(self.original_elem.type, FolderType):
            folderMode = 'true'

        if isinstance(self.original_elem.type, FileStubType):
            folderMode = 'false'

        mainID = self.prefix + "_" + self.get_dest()
        fdElem = etree.SubElement(baseElem, "filedialog", id=mainID, hint=self.get_help())
        fdElem.set('output', outputMode)
        fdElem.set('folder', folderMode)
        #root.set('location', '')
        fdElem.set('multiples', 'true' if self.get_nargs() == '+' else 'false')
        fdElem.set('multiples_delim', ' ')
        fdElem.set('relocateWSL', 'true')
        #root.set('filter', '')
        fdElem.text = self.get_name()


        """
        STARTING WITH EXEC NETWORK HERE
        """
        ifElem = etree.Element("if", id=mainID+"_val")
        ifElem.set('comp', 'IS_SET')
        ifElem.set('value1', "${{{0:}}}".format(mainID))
        ifElem.set('sep', ' ')

        clConst = etree.SubElement(ifElem, 'const')
        clConst.text = self.get_name()

        clValue = etree.SubElement(ifElem, 'value')
        clValue.set('from', "${{{0:}}}".format(mainID))

        #elseElem = etree.SubElement(ifElem, 'else')

        allExecElems = []
        allExecElems.append(ifElem)

        return (ifElem.get('id'), groupElem, allExecElems)


    def make_checkbox(self):

        mainID = self.prefix + "_" + self.get_dest()


        root = etree.Element("checkbox", id=mainID, hint=self.get_help())
        root.set('value', self.get_name())
        root.set('value_unselected', '')
        root.text = self.get_name()

        """
        STARTING WITH EXEC NETWORK HERE
        """
        ifElem = etree.Element("if", id=mainID+"_val")
        ifElem.set('comp', 'IS_SET')
        ifElem.set('value1', "${{{0:}}}".format(mainID))

        clValue = etree.SubElement(ifElem, 'value')
        clValue.set('from', "${{{0:}}}".format(mainID))

        allExecElems = []
        allExecElems.append(ifElem)

        return (ifElem.get('id'), root, allExecElems)


    def make_text_input(self):

        mainID = self.prefix + "_" + self.get_dest()

        groupElem = etree.Element('group')
        baseElem = etree.SubElement(groupElem, "hgroup")

        labelElem = etree.SubElement(baseElem, "label")
        labelElem.text = self.get_dest()


        elemType = "string" #password
        if self.get_type() == int:
            elemType = "int"
        elif self.get_type() == float:
            elemType = "float"

        inputElem = etree.SubElement(baseElem, "input", id=mainID, type=elemType, hint=self.get_help()) # multi=False, min=None, max=None

        defValue = self.get_default()
        if defValue != None:

            if isinstance(defValue, Enum):
                defValue = defValue.name

                inputElem.text = str(defValue)


        """
        STARTING WITH EXEC NETWORK HERE
        """
        ifElem = etree.Element("if", id=mainID+"_val")
        ifElem.set('comp', 'IS_SET')
        ifElem.set('value1', "${{{0:}}}".format(mainID))
        ifElem.set('sep', ' ')


        clConst = etree.SubElement(ifElem, 'const')
        clConst.text = self.get_name()

        clValue = etree.SubElement(ifElem, 'value')
        clValue.set('from', "${{{0:}}}".format(mainID))

        allExecElems = []
        allExecElems.append(ifElem)

        return (ifElem.get('id'), groupElem, allExecElems)

class bioGUIParser:

    @classmethod
    def makeSubparserGroup(cls, subPrefix, parserName):

        root = etree.Element("group", id=subPrefix + "_group")
        root.set('title', parserName)
        root.set('checkable', 'true')
        return root

    @classmethod
    def makeBaseGroup(cls, parserName):

        root = etree.Element("group", id='programgroup')
        root.set('title', parserName)
        root.set('exclusive', 'true')
        return root

    @classmethod
    def addLaunchOption(cls, allParams, executable, progPath, subProg):

        actionButton = etree.Element("action", hint="Click to Run Program", id=subProg + "_run")
        actionButton.set('program', subProg)
        actionButton.text = "Run Program"

        actionNetwork = etree.Element("execute")
        actionNetwork.set('program', subProg)
        actionNetwork.set('exec', executable)
        actionNetwork.set('param', __file__ + " " + subProg + " " + " ".join(allParams))
        actionNetwork.set('wsl', 'true')

        if (progPath):
            actionNetwork.set('location', progPath)

        coutOut = etree.SubElement(actionNetwork, 'output')
        coutOut.set('type', 'COUT')
        coutOut.set('color', 'green')
        coutOut.set('to', 'outstream')

        coutOut = etree.SubElement(actionNetwork, 'output')
        coutOut.set('type', 'CERR')
        coutOut.set('color', 'red')
        coutOut.set('to', 'errstream')

        return (actionButton, [actionNetwork])

    @classmethod
    def makeOutputGroup(cls):

        groupRoot = etree.Element("group")

        streamBox = etree.SubElement(groupRoot, "streambox", id="output1")
        coutStream = etree.SubElement(streamBox,"stream", id="outstream")
        coutStream.text = "COUT"

        cerrStream = etree.SubElement(streamBox,"stream", id="errstream")
        cerrStream.text = "CERR"

        return groupRoot

    @classmethod
    def makeQueryVar(cls, var):
        return  "${{{0:}}}".format(var)


    @classmethod
    def handleSubparser(cls, subCmd, subCmdHelp, subCmdExec, argActions):

        windowParserRoot = cls.makeSubparserGroup(subCmd, subCmd)
        infoLabel = etree.SubElement(windowParserRoot, "label")
        infoLabel.text = subCmdHelp

        parserCLArgs = []
        execElems = []

        for argAction in argActions:

            aelem = ActionElem(argAction, subCmd)
            (aelemID, aelemWindow, aelemExec) = aelem.make_biogui_component()

            windowParserRoot.append(aelemWindow)

            if aelemID != None:
                parserCLArgs.append(cls.makeQueryVar(aelemID))

            if aelemExec != None:
                execElems += aelemExec

        (button, network) = cls.addLaunchOption(parserCLArgs, subCmdExec, None, subCmd)
        windowParserRoot.append(button)

        execElems += network

        return (windowParserRoot, execElems)


    @classmethod
    def makeBioGUIForm(cls, parser):


        templateRoot = etree.Element("template", description="poreSTAT", title="poreSTAT")
        windowRoot = etree.SubElement(templateRoot, "window", title="poreSTAT")
        executionRoot = etree.SubElement(templateRoot, "execution")

        introLabel = etree.SubElement(windowRoot, 'label')
        introLabel.text = parser.format_help()


        allWindowParsers = cls.makeBaseGroup("poreStat")
        allExecutionParsers = None

        for spAction in parser._subparsers._group_actions:

            if spAction.choices == None:

                pass
            else:

                for (subcmd, choiceInputs) in spAction.choices.items():

                    subActions = choiceInputs._actions
                    (windowParserRoot, network) = cls.handleSubparser(subcmd, choiceInputs.format_help(), 'python3', subActions)

                    for x in network:
                        executionRoot.append(x)

                    allWindowParsers.append(windowParserRoot)


        windowRoot.append(allWindowParsers)

        if allExecutionParsers != None:
            executionRoot.append(allExecutionParsers)

        outGroup = cls.makeOutputGroup()
        windowRoot.append(outGroup)

        print(etree.tostring(templateRoot, pretty_print=True).decode())

    @classmethod
    def makeBioGUIFormFromTool(cls, tool):

        templateRoot = etree.Element("template", description="poreSTAT", title="poreSTAT")
        windowRoot = etree.SubElement(templateRoot, "window", title="poreSTAT")
        executionRoot = etree.SubElement(templateRoot, "execution")

        introLabel = etree.SubElement(windowRoot, 'label')
        introLabel.text = parser.format_help()


        allWindowParsers = cls.makeBaseGroup("poreStat")
        allExecutionParsers = None

        allInputArgs = []
        for arg in tool.inputs.values():
            arg.prefix = args.prefix
            allInputArgs.append(Argument(arg, 'r'))

        allOutputArgs = []
        for arg in tool.outputs.values():
            arg.prefix = args.prefix
            allOutputArgs.append(Argument(arg, 'w'))


        (windowParserRoot, execElems) = cls.handleSubparser("", "", "", allInputArgs+allOutputArgs)
        allWindowParsers.append(windowParserRoot)

        windowRoot.append(allWindowParsers)

        if execElems != None:
            for x in execElems:
                executionRoot.append(x)

        outGroup = cls.makeOutputGroup()
        windowRoot.append(outGroup)

        print(etree.tostring(templateRoot, pretty_print=True).decode())




if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--cwl', nargs='+', type=argparse.FileType('r'), help='CWL file to read in', required=True)
    parser.add_argument('-o', '--output', type=FolderType('w'), help="Location to store biogui templates in", required=False)
    parser.add_argument('-p', '--prefix', type=str, help="prefix for all variables", required=False, default="")

    args = parser.parse_args()

    for file in args.cwl:

        argument_names = []

        try:
            tool = Tool(file.name)
        except yaml.scanner.ScannerError:
            sys.exit('File {0} is corrupted or not a CWL tool definition')

        bioGUIParser.makeBioGUIFormFromTool(tool)

        filename = os.path.basename(file.name)
        (base,ext) = os.path.splitext(filename)

        print(base)
