import itertools


# List of lines
def line_cycle():
    return itertools.cycle(("-", "--", ":", "-."))


# List of markers
def marker_cycle():
    return itertools.cycle(('o', 'v', 's', 'p', 'X', '^', '<', '>'))


# List of colors
def color_cycle():
    return itertools.cycle(('C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7'))
    # return itertools.cycle(('b', 'g', 'r', 'c', 'm', 'y', 'k', 'w'))
