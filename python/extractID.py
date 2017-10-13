import sys
from collections import defaultdict

ref = sys.argv[1]
query = sys.argv[2]

refDict = defaultdict(list)

for line in open(ref):
    if line[0] == ">":
        key = line.split(' ')[0][1:]
    else:
        refDict[key].append(line.strip())

for key, value in refDict.items():
    refDict[key] = ''.join(value)

for each in open(query):
    result = refDict[each.strip()]
    print(">" + each + str(result))

