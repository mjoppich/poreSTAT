import sys
from itertools import chain
from collections import defaultdict, Counter


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)



def mergeDicts( dict1, dict2):
    dict3 = {}
    for k, v in chain(dict1.items(), dict2.items()):

        if k in dict3:

            if not type(v) == type(dict3[k]):
                raise Exception("You try to merge two different objects!")

            if type(v) == list:

                dict3[k] = dict3[k] + v

            elif type(v) == set:

                dict3[k] = dict3[k].union(v)

            elif type(v) == dict:

                dict3[k] = mergeDicts(dict3[k], v)

            else:

                dict3[k] = dict3[k] + v
        else:

            dict3[k] = v

    return dict3

def mergeCounter(self, counter1, counter2):

    mergedCounter = Counter()

    for x in counter1:
        mergedCounter[x] = counter1[x]

    for x in counter2:
        mergedCounter[x] += counter2[x]

    return mergedCounter