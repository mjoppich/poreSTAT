import os

class PTToolInterface:

    def __init__(self, parser, subparser):

        self.parser = parser
        self.subparser = subparser


    def print_usage(self):
        self.subparser.print_usage()

    def print_help(self):
        self.subparser.print_help()

    def exec(self, args):
        pass

    def _getsubdirs(self, path):

        alldirs = [o for o in os.listdir(path) if os.path.isdir(os.path.join(path, o))]

        return alldirs


    def manage_folders_reads(self, args):

        if (args.folders == None and args.reads == None):
            print("error: Either folders or reads must be set!")
            self.subparser.print_help()
            exit(-1)

        if args.folders != None:
            folders = args.folders
        else:
            folders = []

        if args.reads != None:

            # if downloads folder, then take the downloads folder + fail/pass plus any batch folder within
            downloadpath = os.path.join(args.reads, "downloads")
            if os.path.isdir( downloadpath ):

                for x in ['fail', 'skip', 'pass']:

                    if os.path.isdir(os.path.join(downloadpath, x)):

                        curfolder = os.path.join(downloadpath, x)

                        folders.append( curfolder )

                        all_subdirs = self._getsubdirs(curfolder)

                        for x in all_subdirs:

                            if x.startswith("batch"):
                                folders.append( os.path.join(curfolder, x) )


            else:
                # if not downloads folder
                folders.append( os.path.join(args.reads, 'fail') )
                folders.append(os.path.join(args.reads, 'skip'))
                folders.append(os.path.join(args.reads, 'pass'))

        return folders