from matplotlib import pyplot as plt


filename = "delay_60_5.txt"

with open(filename) as f:
    filecontent = f.read()
filecontent = filecontent.split("\n")
keys = []
values = []
for item in filecontent:
    if len(item) > 0:
        pair = item.split(" ")
        keys.append(int(pair[0]))
        values.append(int(pair[1]))
values.sort()
fig = plt.figure()
A = range(len(values)+1)
A = A[1:]
plt.plot(A,values)
plt.show()
