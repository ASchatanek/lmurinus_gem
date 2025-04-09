import collections as co

cnt = co.Counter()

thislist = ["w", "blue", "blue", "blue"]

for word in thislist:

    cnt[word] += 1

cnt

for x in cnt.values():

    print(x)
