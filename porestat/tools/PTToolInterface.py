
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