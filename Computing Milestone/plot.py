import matplotlib.pyplot as plt

data = open("RK4.txt", "r")
string = data.read()
data.close()

split = string.split("\n")
x= []
y = []

steps = 200000000 * 50
skips = steps / 5000
dt = 0.0005

for i in range(len(split)):
    if len(split[i]) > 0:
        info = split[i].split(", ")
        y.append(float(info[0]))
        x.append((dt * (i * skips + 1)) / (60 * 60 * 24))
        

plt.plot(x, y)
plt.xlim(xmin=0)
plt.ylim(ymin = 0)
plt.xlabel("Time, days")
plt.ylabel("Distance From L2, m")
plt.title("Runge-Kutta")
plt.show()
