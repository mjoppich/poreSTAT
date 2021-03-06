import os

from porestat.utils import eprint


class ArgObj(object):
    pass

from abc import ABC, abstractmethod


class PSToolInterfaceFactory(ABC):

    def __init__(self, parser, subparser, which):

        self.parser = parser
        self.subparser = subparser
        self.which = which

        self.subparser.set_defaults(which=which)


    def print_usage(self):
        self.subparser.print_usage()

    def print_help(self):
        self.subparser.print_help()

    def _makeArguments(self, args):

        newArgs = ArgObj()

        for x in vars(args):
            newArgs.__dict__[x] = args.__dict__[x]

        return newArgs


    @abstractmethod
    def _addParser(self, subparsers, which):
        pass


class PSToolException(Exception):

    def __init__(self, msg):

        super(PSToolException, self).__init__()

        self.msg = msg

    def __str__(self):
        fp = super().__str__()
        fp += "\n\n"
        fp += self.msg

        return fp

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

    @classmethod
    def hasArgument( self, argName, args ):
        return argName in args.__dict__

    def makeKey(self, runUserName, args, runid):
        return ",".join(runUserName) if self.hasArgument('user_run', args) and args.user_run else runid

    def manage_folders_reads(self, args):

        if (args.folders == None and args.reads == None and args.mreads == None):
            eprint("error: Either folders or reads must be set!")
            raise PSToolException("error: Either folders or reads must be set!")

        if args.folders != None:
            folders = args.folders
        else:
            folders = []

        if args.mreads != None:

            for mread in args.mreads:
                folders.append(mread)

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