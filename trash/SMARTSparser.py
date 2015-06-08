import re

s = 'C1=CC2=C(C=C1)C=CC=C2'


def brainfuck2list(brainfuck):
  while brainfuck:               #if list is empty then finish
    e = brainfuck.pop(0)
    if e not in ("(", ")"):
      yield e
    elif e == "(":
      yield list(brainfuck2list(brainfuck))
    else:
      break

print([_ for _ in brainfuck2list([x for x in re.split('(\(|\))', s) if x != ''])])
