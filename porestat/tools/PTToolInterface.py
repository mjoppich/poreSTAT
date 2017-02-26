
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
            folders.append( '' )

        return folders