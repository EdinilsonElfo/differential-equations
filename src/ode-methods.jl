const euler = ExplicitEuler()
const midp = MidPoint()
const heun = Heun()

a = [
    0.0 0.0
    0.5 0.0
]

b = [0.0, 1.0]

c = [0.0, 0.5]

const rk2 = RungeKutta(a, b, c, 2)


a = [
    0.0 0.0 0.0
    1/3 0.0 0.0
    0.0 2/3 0.0
]

b = [1/4, 0.0, 3/4]

c = [0.0, 1/3, 2/3]

const rk3 = RungeKutta(a, b, c, 3)

a = [
    0.0 0.0 0.0 0.0
    0.5 0.0 0.0 0.0
    0.0 0.5 0.0 0.0
    0.0 0.0 1.0 0.0
]

b = [1/6, 1/3, 1/3, 1/6]

c = [0.0, 0.5, 0.5, 1.0]

const rk4 = RungeKutta(a, b, c, 4)

