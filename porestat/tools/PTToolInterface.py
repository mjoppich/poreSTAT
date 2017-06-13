import os

from porestat.utils import eprint


class ArgObj(object):
    pass

class PSToolInterfaceFactory:

    def __init__(self, parser, subparser):

        self.parser = parser
        self.subparser = subparser


    def print_usage(self):
        self.subparser.print_usage()

    def print_help(self):
        self.subparser.print_help()

    def _makeArguments(self, args):

        newArgs = ArgObj()

        for x in vars(args):
            newArgs.__dict__[x] = args.__dict__[x]

        return newArgs



class PSToolException(Exception):

    def __init__(self, msg):

        super(PSToolException, self).__init__()

        self.msg = msg

class PSToolInterface:

    def __init__(self, args):

        self.args = args

    def exec(self, args):
        pass

    def _getsubdirs(self, path):

        alldirs = [o for o in os.listdir(path) if os.path.isdir(os.path.join(path, o))]

        return alldirs

    def _makeArguments(self, args):

        newArgs = ArgObj()

        for x in vars(args):
            newArgs.__dict__[x] = args.__dict__[x]

        return newArgs

    def hasArgument( self, argName, args ):

        return argName in args.__dict__

    def manage_folders_reads(self, args):

        if (args.folders == None and args.reads == None):
            eprint("error: Either folders or reads must be set!")
            raise PSToolException("bla")

        if args.folders != None:
            folders = args.folders
        else:
            folders = []

        if args.reads != None:

            for readfolder in args.reads:

                # if downloads folder, then take the downloads folder + fail/pass plus any batch folder within
                downloadpath = os.path.join(readfolder, "downloads")
                if os.path.isdir(downloadpath):

                    for x in ['fail', 'skip', 'pass']:

                        if os.path.isdir(os.path.join(downloadpath, x)):

                            curfolder = os.path.join(downloadpath, x)

                            folders.append(curfolder)

                            all_subdirs = self._getsubdirs(curfolder)

                            for x in all_subdirs:

                                if x.startswith("batch"):
                                    folders.append(os.path.join(curfolder, x))


                else:

                    for x in os.walk(readfolder):
                        eprint("Adding path: " + x[0])

                        folders.append(x[0])

        return folders