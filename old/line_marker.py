import itertools

# List of lines
lines = itertools.cycle(("-", "--", ":", "-."))

# List of markers
markers = itertools.cycle(('o', 'v', 's', 'p', 'X', '^', '<', '>'))

# List of colors
# colors = itertools.cycle(('b', 'g', 'r', 'c', 'm', 'y', 'k', 'w'))
colors = itertools.cycle(('C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7'))
