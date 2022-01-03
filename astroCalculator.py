import math

def form_1():
    c = 299792458
    print("Please enter the values")
    lambdaVal = input("\u03BB's value: ")
    f = input("frequency value: ")
    if f == "?":
        frequency = c/float(lambdaVal)
        f = float(frequency)
        print("The answer is f = " + str(f) + "\n")
    else:
        lamb = c/float(f)
        lambdaVal = float(lamb)
        print("The answer is \u03BB = " + str(lambdaVal) + "\n")


def form_2():
    h = 6.626 * math.pow(10, -34)
    c = 299792458
    print("You can choose one of the formulas:\n"
          "1- E=h*c/\u03BB\n"
          "2- E=hf\n")
    formula_choice = input("Please write the number of the formula that you chose: ")
    if formula_choice =="1":
        print("Please enter the values")
        lambdaVal = input("\u03BB's value: ")
        evalue = input("E's value: ")
        if evalue == "?":
            energy = h*c/float(lambdaVal)
            evalue = float(energy)
            print("The answer is E = " + str(evalue) + "\n")
        else:
            lamb = h*c/float(evalue)
            lambdaVal = float(lamb)
            print("The answer is \u03BB = " + str(lambdaVal) + "\n")
    else:
        print("Please enter the values")
        f = input("f's value: ")
        evalue = input("E's value: ")
        if evalue == "?":
            energy = h * float(f)
            evalue = float(energy)
            print("The answer is E = " + str(evalue) + "\n")
        else:
            freq = float(evalue)/h
            f = float(freq)
            print("The answer is f = " + str(f) + "\n")

def form_3():
    c = 299792458
    print("Please enter the values:\n")
    deltLamb = input("\u0394\u03BB's value: ")
    lamb  =input("\u03BB's value: ")
    speed  = input("v's value: ")
    if deltLamb == "?":
        delta = float(lamb)*(float(speed)/c)
        deltLamb = float(delta)
        print("The answer is \u0394\u03BB = " + str(deltLamb) + "\n")
    elif lamb == "?":
        lambd = (float(deltLamb)*c)/float(speed)
        lamb = float(lambd)
        print("The answer is \u03BB = " + str(lamb) + "\n")
    else:
        v = (float(deltLamb)*c)/float(speed)
        speed = float(v)
        print("The answer is v = " + str(speed) + "\n")

def form_4():
    G = 6.67 * math.pow(10, -11)
    print("Please enter the values\n")
    force = input("F's value: ")
    bigM = input("M's value: ")
    m = input("m's value: ")
    r = input("r's value: ")
    if force == "?":
        f = (G*float(bigM)*float(m))/math.pow(float(r),2)
        force = float(f)
        print("The answer is F = " + str(force) + "\n")
    elif bigM == "?":
        massM = (float(force)* math.pow(float(r),2))/(float(m)*G)
        bigM = float(massM)
        print("The answer is M = " + str(bigM) + "\n")
    elif m == "?":
        smallM = (float(force)* math.pow(float(r),2))/(float(bigM)*G)
        m = float(smallM)
        print("The answer is m = " + str(m) + "\n")
    else:
        rad = math.sqrt((G*float(bigM)*float(m))/float(force))
        r = float(rad)
        print("The answer is r = " + str(r) + "\n")

def form_5():
    G = 6.67 * math.pow(10, -11)
    print("Please enter the values\n")
    bigM = input("M's value: ")
    a = input("a's value: ")
    r = input("r's value: ")
    if a == "?":
        accel = (G*float(bigM))/math.pow(float(r),2)
        a = float(accel)
        print("The answer is a = " + str(a) + "\n")
    elif bigM == "?":
        massM = (float(a)* math.pow(float(r),2))/(G)
        bigM = float(massM)
        print("The answer is M = " + str(bigM) + "\n")
    else:
        rad = math.sqrt((G * float(bigM)) / float(a))
        r = float(rad)
        print("The answer is r = " + str(r) + "\n")

def form_6():
    G = 6.67 * math.pow(10, -11)
    print("Please enter the values\n")
    bigM = input("M's value: ")
    v = input("V_e's value: ")
    r = input("r's value: ")
    if v == "?":
        vel = math.sqrt((2*G*float(bigM))/float(r))
        v = float(vel)
        print("The answer is V_e = " + str(v) + "\n")
    elif bigM == "?":
        massM = (float(r)* math.pow(float(v),2))/(2*G)
        bigM = float(massM)
        print("The answer is M = " + str(bigM) + "\n")
    else:
        rad = (2*G * float(bigM)) / math.pow(float(v),2)
        r = float(rad)
        print("The answer is r = " + str(r) + "\n")

def form_7():
    G = 6.67 * math.pow(10, -11)
    c = 299792458
    print("You can choose one of the formulas:\n"
          "1- R = 2GM/c^2\n"
          "2- R = 3*(M/M_sun)\n")
    formula_choice = input("Please write the number of the formula that you chose: ")
    if formula_choice =="1":
        print("Please enter the values")
        bigM = input("M's value: ")
        r = input("R's value: ")
        if r == "?":
            rad = (2*G*float(bigM))/math.pow(c,2)
            r = float(rad)
            print("The answer is R_EH = " + str(r) + "\n")
        else:
            greatM = (float(r)*math.pow(c,2))/(2*G)
            bigM = float(greatM)
            print("The answer is M = " + str(bigM) + "\n")
    else:
        print("Please enter the values")
        m = input("M's value: ")
        mSun = input("M_sun's value: ")
        r = input("R's value: ")
        if r == "?":
            rad = 3*(float(m)/float(mSun))
            r = float(rad)
            print("The answer is R_s = " + str(r) + "\n")
        elif m == "?":
            bigM = (float(mSun)*float(r))/3
            m = bigM
            print("The answer is M = " + str(m) + "\n")
        else:
            bigMsun = (3*float(m))/float(r)
            mSun = bigMsun
            print("The answer is M_sun = " + str(mSun) + "\n")

def form_8():
    print("These peak wavelengths are directly related to the surface temperature of a star.\n")
    print("Please enter the values\n")
    temp = input("T's value: ")
    lambdaVal = input("\u03BB's value: ")

    if lambdaVal == "?":
        lamVal = 0.0029/float(temp)
        lambdaVal = float(lamVal)
        print("The answer is \u03BB_peak = " + str(lambdaVal) + "\n")
    else:
        tempt = 0.0029/float(lambdaVal)
        temp = float(tempt)
        print("The answer is T = " + str(temp) + "\n")

def form_9():
    c = 299792458
    print("The mass is not really lost, it has transformed into thermal energy, as described by Einstein's famous equation.\n")
    print("Please enter the values\n")
    energy = input("E's value: ")
    mass = input("m's value: ")

    if energy == "?":
        eng = float(mass)*math.pow(c,2)
        energy = float(eng)
        print("The answer is E = " + str(energy) + "\n")
    else:
        m = float(energy)/math.pow(c,2)
        mass = float(m)
        print("The answer is m = " + str(mass) + "\n")

def form_10():
    c = 299792458
    print("Any nucleus is made up of N-numbers neutrons and Z-numbers of protons where N and Z depend on the element. The binding energy is dened by adding the mass of all the protons and all the neutrons and subtracting the mass of the nucleus and multiplying by c-squared.\n")
    print("Please enter the values\n")
    energy = input("E's value: ")
    z = input("Z's value: ")
    n = input("N's value: ")
    m = input("M_nucleus's value: ")
    if energy == "?":
        eng = (float(z)+float(n)-float(m))*math.pow(c,2)
        energy = float(eng)
        print("The answer is E = " + str(energy) + "\n")
    elif z == "?":
        z_mp = (float(energy)/math.pow(c,2))-float(n)+float(m)
        z = float(z_mp)
        print("The answer is Z = " + str(z) + "\n")
    elif n == "?":
        n_mp = (float(energy)/math.pow(c,2))-float(z)+float(m)
        n = float(n_mp)
        print("The answer is N = " + str(n) + "\n")
    else:
        m_nuc = -(float(energy) / math.pow(c, 2)) + float(z) + float(n)
        m = float(m_nuc)
        print("The answer is M = " + str(m) + "\n")


def form_11():
    c = 299792458
    print("Please enter the values\n")
    haveT0 = input("Do we have the values of t_0?(y/n) ")
    if haveT0 == "n":
        d = input("d's value: ")
        v = input("v's value: ")
        tZero = (2*float(d))/c
        t = float(tZero)/math.sqrt(1-math.pow((float(v)/c),2))
        print("The answer is t_0 = " + str(tZero) + "\n")
        print("The answer is t = " + str(t) + "\n")
    else:
        d = input("d's value: ")
        v = input("v's value: ")
        tZero = input("t_0's value: ")
        t = input("t's value: ")

        if t == "?":
            t = float(tZero)/math.sqrt(1-math.pow((float(v)/c),2))
            print("The answer is t = " + str(t) + "\n")
        elif d == "?":
            dist = float(tZero)*c/2
            print("The answer is d = " + str(d) + "\n")
        else:
            speed = c*math.sqrt(1-(math.pow((float(tZero)/float(t)),2)))
            v = float(speed)
            print("The answer is v = " + str(v) + "\n")


def form_12():
    c = 299792458
    print("Please enter the values\n")
    v = input("v's value: ")
    lZero = input("L_0's value: ")
    l = input("L's value: ")
    if l == "?":
        length = float(lZero)*math.sqrt(1-math.pow((float(v)/c),2))
        l = float(length)
        print("The answer is L = " + str(l) + "\n")
    elif lZero == "?":
        lenZero = float(l)/math.sqrt(1-math.pow((float(v)/c),2))
        lZero = float(lenZero)
        print("The answer is L_0 = " + str(lZero) + "\n")
    else:
        speed = c*math.sqrt(1-math.pow((float(l)/float(lZero)),2))
        v = float(speed)
        print("The answer is v = " + str(v) + "\n")

def form_13():
    c = 299792458
    print("Please enter the values\n")
    gamma = input("\u03B3's value: ")
    v = input("v's value: ")
    if gamma == "?":
        loren = 1/(math.sqrt(1-math.pow((float(v)/c),2)))
        gamma = float(loren)
        print("The answer is \u03B3 = " + str(gamma) + "\n")
    else:
        speed = c*math.sqrt(1-math.pow((1/float(gamma)),2))
        v = float(speed)
        print("The answer is v = " + str(v) + "\n")


def form_14():
    c = 299792458
    G = 6.67 * math.pow(10, -11)
    print("Please enter the values\n")
    tPlanet = input("\u0394t_planet's value: ")
    tOrbit = input("\u0394t_orbit's value: ")
    mass = input("M's value: ")
    r = input("R's value: ")
    if tPlanet == "?":
        planet = float(tOrbit)*(math.sqrt(1-((2*G*float(mass))/(float(r)*math.pow(c,2)))))
        tPlanet = float(planet)
        print("The answer is \u0394t_planet = " + str(tPlanet) + "\n")
    elif tOrbit == "?":
        orbit = float(tPlanet)/(math.sqrt(1 - ((2 * G * float(mass)) / (float(r) * math.pow(c, 2)))))
        tOrbit = float(orbit)
        print("The answer is \u0394t_orbit = " + str(tOrbit) + "\n")
    elif mass == "?":
        bigM = ((float(r)*math.pow(c,2))/2*G) * (1-math.pow((float(tPlanet)/float(tOrbit)),2))
        mass =float(bigM)
        print("The answer is M = " + str(mass) + "\n")
    else:
        rad = (2*G*float(mass))/(math.pow(c,2)*(1-math.pow((float(tPlanet)/float(tOrbit)),2)))
        r = float(rad)
        print("The answer is R = " + str(r) + "\n")


def form_15():
    print("Please enter the values\n")
    separate = input("Do we have the values of M_1 and M_2 or one of them?(y/n)\n")
    if separate == "y":
        m1 = input("M_1's value: ")
        m2 = input("M_2's value: ")
        a = input("a's value: ")
        p = input("p's value: ")
        if m1 == "?":
            mOne = (math.pow(float(a),3)/math.pow(float(p),2))-float(m2)
            m1 = float(mOne)
            print("The answer is M_1 = " + str(m1) + "\n")
        elif m2 == "?":
            mTwo = (math.pow(float(a), 3) / math.pow(float(p), 2)) - float(m1)
            m2 = float(mTwo)
            print("The answer is M_2 = " + str(m2) + "\n")
        elif a == "?":
            dist = math.pow((float(m1)+float(m2))/math.pow(float(p),2),1/3)
            a = float(dist)
            print("The answer is a = " + str(a) + "\n")
        else:
            period = math.sqrt(math.pow(float(a), 3)/(float(m1) + float(m2)))
            p = float(period)
            print("The answer is p = " + str(p) + "\n")
    else:
        bigM = input("value of M_1+M_2: ")
        a = input("a's value: ")
        p = input("p's value: ")
        if bigM == "?":
            mSum = (math.pow(float(a), 3) / math.pow(float(p), 2))
            bigM = float(mSum)
            print("The answer is M_1 + M_2 = " + str(bigM) + "\n")
        elif a == "?":
            dist = math.pow(float(bigM)/math.pow(float(p),2),1/3)
            a = float(dist)
            print("The answer is a = " + str(a) + "\n")
        else:
            period = math.sqrt(math.pow(float(a), 3)/float(bigM))
            p = float(period)
            print("The answer is p = " + str(p) + "\n")

def form_16():
    print("Please enter the values\n")
    p = input("p's value: ")
    m = input("m's value: ")
    v = input("v's value: ")
    if p == "?":
        momentum = float(m)*float(v)
        p = float(momentum)
        print("The answer is p = " + str(p) + "\n")
    elif m == "?":
        mass = float(p) / float(v)
        m = float(mass)
        print("The answer is a = " + str(m) + "\n")
    else:
        speed = float(p) / float(m)
        v = float(speed)
        print("The answer is v = " + str(v) + "\n")

def form_17():
    print("Please enter the values\n")
    l = input("L's value: ")
    i = input("I's value: ")
    omega = input("\u03C9's value: ")
    if l == "?":
        linear = float(i)*float(omega)
        l = float(linear)
        print("The answer is L = " + str(l) + "\n")
    elif i == "?":
        inertia = float(l) / float(omega)
        i = float(inertia)
        print("The answer is I = " + str(i) + "\n")
    else:
        omg = float(l) / float(i)
        omega = float(omg)
        print("The answer is \u03C9 = " + str(omega) + "\n")

def form_18():
    G = 6.67 * math.pow(10, -11)
    print("Please enter the values\n")
    ep = input("E's value: ")
    bigM = input("M's value: ")
    tinniM = input("m's value: ")
    r = input("R's value: ")
    if ep == "?":
        potential = (G*float(bigM)*float(tinniM))/float(r)
        ep = float(potential)
        print("The answer is E_potential = " + str(ep) + "\n")
    elif r == "?":
        dist = (G*float(bigM)*float(tinniM))/float(ep)
        r = float(dist)
        print("The answer is R = " + str(r) + "\n")
    elif bigM == "?":
        bigMass = (float(ep)*float(r)) / float(tinniM)*G
        bigM = float(bigMass)
        print("The answer is M = " + str(bigM) + "\n")
    else:
        minMass = (float(ep) * float(r)) / float(bigM) * G
        tinniM = float(minMass)
        print("The answer is m = " + str(tinniM) + "\n")

def form_19():
    G = 6.67 * math.pow(10, -11)
    c = 299792458
    pi = 3.14
    print("Please enter the values\n")
    edd = input("L's value: ")
    bigM = input("M's value: ")
    sigm = input("\u03C3's value: ")
    m = input("m_p's value: ")
    if edd == "?":
        lEdd = (4*pi*G*float(bigM)*float(m)*c)/float(sigm)
        edd = float(lEdd)
        print("The answer is L_edd = " + str(edd) + "\n")
    elif sigm == "?":
        sigmas = (4*pi*G*float(bigM)*float(m)*c)/float(edd)
        sigm = float(sigmas)
        print("The answer is \u03C3 = " + str(sigm) + "\n")
    elif bigM == "?":
        bigMass = float(edd)*float(sigm)/ (4 * pi * G * float(m) * c)
        bigM = float(bigMass)
        print("The answer is M = " + str(bigM) + "\n")
    else:
        minMass = float(edd)*float(sigm)/ (4 * pi * G * float(bigM) * c)
        tinniM = float(minMass)
        print("The answer is m = " + str(tinniM) + "\n")

def form_20():
    G = 6.67 * math.pow(10, -11)
    print("Please enter the values\n")
    a = input("a's value: ")
    bigM = input("M's value: ")
    hight = input("h's value: ")
    r = input("r's value: ")
    if a == "?":
        accel = (2*G*float(bigM)*float(hight))/math.pow(float(r),3)
        a = float(accel)
        print("The answer is a = " + str(a) + "\n")
    elif r == "?":
        dist = math.pow((2*G*float(bigM)*float(hight))/float(a),1/3)
        r = float(dist)
        print("The answer is R = " + str(r) + "\n")
    elif bigM == "?":
        bigMass = math.pow(float(r),3)*float(a)/(2*G*float(hight))
        bigM = float(bigMass)
        print("The answer is M = " + str(bigM) + "\n")
    else:
        h = math.pow(float(r),3)*float(a)/(2*G*float(bigM))
        hight = float(h)
        print("The answer is h = " + str(hight) + "\n")


def form_21():
    print("Please enter the values\n")
    t = input("T_fall's value: ")
    bigM = input("M's value: ")
    mSun = input("M_sun's value: ")

    if t == "?":
        fall = 15*math.pow(10,-6)*(float(bigM)/float(mSun))
        t = float(fall)
        print("The answer is T_fall = " + str(t) + "\n")
    elif bigM == "?":
        bigMass = float(t)*float(mSun)/ 15*math.pow(10,-6)
        bigM = float(bigMass)
        print("The answer is M = " + str(bigM) + "\n")
    else:
        sun = 15*math.pow(10,-6)*(float(bigM)/float(t))
        mSun = float(sun)
        print("The answer is M_sun = " + str(mSun) + "\n")

def form_22():
    G = 6.67 * math.pow(10, -11)
    c = 299792458
    print("Please enter the values\n")
    a = input("a's value: ")
    bigM = input("M's value: ")
    j = input("J's value: ")

    if a == "?":
        rotate = (float(j)*c)/(math.pow(float(bigM),2)*G)
        a = float(rotate)
        print("The answer is a = " + str(a) + "\n")
    elif bigM == "?":
        bigMass = math.sqrt((float(j)*c)/(float(a)*float(j)))
        bigM = float(bigMass)
        print("The answer is M = " + str(bigM) + "\n")
    else:
        angular = (float(a)*math.pow(float(bigM),2)*G)/c
        j = float(angular)
        print("The answer is J = " + str(j) + "\n")

def form_23():

    print("Please enter the values\n")
    f = input("F_max's value: ")
    bigM = input("M's value: ")
    mSun = input("M_sun's value: ")
    if f == "?":
        frequency = 16000 * (float(mSun)/float(bigM))
        f = float(frequency)
        print("The answer is F_max = " + str(f) + "\n")
    elif bigM == "?":
        bigMass = 16000 * (float(mSun)/float(f))
        bigM = float(bigMass)
        print("The answer is M = " + str(bigM) + "\n")
    else:
        sun = (float(bigM)*float(f))/16000
        mSun = float(sun)
        print("The answer is M_sun = " + str(mSun) + "\n")

def form_24():
    h = 6.6*math.pow(10,-34)
    print("Please enter the values\n")
    lamb = input("\u03BB's value: ")
    v = input("v's value: ")
    m = input("m's value: ")
    if lamb == "?":
        lambd = h/(float(m)*float(v))
        lamb = float(lambd)
        print("The answer is \u03BB = " + str(lamb) + "\n")
    elif m == "?":
        mass = h/(float(lamb)*float(v))
        m = float(mass)
        print("The answer is m = " + str(m) + "\n")
    else:
        speed = h/(float(m)*float(lamb))
        v = float(speed)
        print("The answer is v = " + str(v) + "\n")

def form_25():
    k = 1.38064852*math.pow(10,-23)
    print("Please enter the values\n")
    s = input("S's value: ")
    w = input("W's value: ")
    if s == "?":
        entropy = k*math.log10(float(w))
        s = float(entropy)
        print("The answer is S = " + str(s) + "\n")
    else:
        equip = math.pow(10,(float(s)/k))
        w = float(equip)
        print("The answer is W = " + str(w) + "\n")

def form_26():
    k = 1.38064852*math.pow(10,-23)
    l = 1.6*math.pow(10,-35)
    print("Please enter the values\n")
    s = input("S's value: ")
    a = input("A's value: ")
    if s == "?":
        entropy = (k*float(a))/(4*math.pow(l,2))
        s = float(entropy)
        print("The answer is S_BH = " + str(s) + "\n")
    else:
        area = (4*float(s)*math.pow(l,2))/k
        a = float(area)
        print("The answer is A = " + str(a) + "\n")

def form_27():
    print("Please enter the values\n")
    t = input("T's value: ")
    k = input("k's value: ")

    if t == "?":
        temp = float(k)/(2*3.14)
        t = float(temp)
        print("The answer is T = " + str(t) + "\n")
    else:
        kappa = 2*3.14*float(t)
        k = float(kappa)
        print("The answer is k = " + str(k) + "\n")

def form_28():
    G = 6.67 * math.pow(10, -11)
    c = 299792458
    pi = 3.14
    print("Please enter the values\n")
    t = input("T's value: ")
    k = input("k's value: ")
    m = input("M's value: ")
    hBar = input("h_bar's value: ")
    if t == "?":
        temp = (float(hBar)*math.pow(c,3))/(8*pi*G*float(m)*float(k))
        t = float(temp)
        print("The answer is T = " + str(t) + "\n")
    elif k == "?":
        kappa = (float(hBar)*math.pow(c,3))/(8*pi*G*float(m)*float(t))
        k = float(kappa)
        print("The answer is k = " + str(k) + "\n")
    elif m == "?":
        mass = (float(hBar) * math.pow(c, 3)) / (8 * pi * G * float(t) * float(t))
        m = float(mass)
        print("The answer is M = " + str(m) + "\n")
    else:
        h = (float(t)*8*G*pi*float(m)*float(k))/math.pow(c,3)
        hBar = float(h)
        print("The answer is h_Bar = " + str(h) + "\n")


def form_29():
    print("Please enter the values\n")
    t = input("T's value: ")
    lamb = input("\u03BB_peak's value: ")
    w = input("W's value: ")

    if t == "?":
        temp = float(w)/float(lamb)
        t = float(temp)
        print("The answer is T = " + str(t) + "\n")
    elif w == "?":
        waves = float(t)*float(lamb)
        w = float(waves)
        print("The answer is W = " + str(w) + "\n")
    else:
        lambd = float(w)/float(t)
        lamb = float(lambd)
        print("The answer is \u03BB_peak = " + str(lamb) + "\n")

def form_30():
    G = 6.67 * math.pow(10, -11)
    c = 299792458
    print("Please enter the values\n")
    bigM = input("M's value: ")
    r = input("R's value: ")

    if r == "?":
        shadow = 2.6*((2*G*float(bigM))/math.pow(c,2))
        r = float(shadow)
        print("The answer is R_shadow = " + str(r) + "\n")
    else:
        mass = (float(r)*math.pow(c,2))/(2.6*2*G)
        bigM = float(mass)
        print("The answer is M = " + str(bigM) + "\n")
def main():

    key = True
    while key:
        print("*****Welcome to the Astro Calculator!*****\n")
        print("Please choose the number of the formula that you would like to use. Please put '?' in front of the variable that you are looking for when you are using a "
              "formula.\n")

        print(""
              "0- Exit\n"
              "1- The equation relating frequency and wavelength of a photon to the speed of light.(\u03BB*f=c)\n"
              "2- The photon energy.(E=h*c/\u03BB or E=hf)\n"
              "3- Doppler Shift Formula.(\u0394\u03BB=\u03BB*v/c)\n"
              "4- Newton's Universal Law Of Gravitation.(F=GMm/r^2)\n"
              "5- Formula for calculating the acceleration due to gravity. (a=GM/r^2)\n"
              "6- Escape Velocity.(v = sqrt(2GM/r))\n"
              "7- Schwarzschild Radius Formulas.(R=2GM/c^2 or R = 3*(M/M_sun)\n"
              "8- Wien's Law.(\u03BB=0.0029/T)\n"
              "9- Einstein's Famous Equation.(E=mc^2)\n"
              "10- Binding Energy.(E = Z_mp+N_mp-m_nucleus)c^2)\n"
              "11- Time Dilation.(t_0 = 2d/c, t = sqrt(t_0/(1-(v/c)^2)))\n"
              "12- Length Contraction.(L=L_0*sqrt(1-(v/c)^2)\n"
              "13- Lorentz Factor or Gamma Factor.(\u03B3(v)=1/sqrt(1-(v/c)^2))\n"
              "14- Geodesics. To calculate how time is wrapped in a strong gravitational field.(\u0394t_planet = \u0394t_orbit * sqrt(1-(2GM/Rc^2)))\n"
              "15- Kepler's third law is an equation that relates the masses of the stars to the orbital period and the total distance between the stars.(M_1+M_2=a^3/p^2)\n"
              "16- Momentum.(p = m * v)\n"
              "17- That way, angular momentum, L, can be expressed in a similar way to linear momentum as the product of Moment of Inertia, I, multiplied by the Angular velocity, omega.(L = I * \u03C9)\n"
              "18- A black hole with mass M and radius R, and a test particle with a mass of little-m, the gravitational potentialenergy for the infalling particle is: E = GMm/R\n"
              "19- Eddington limit. The Eddington limit describes a natural limit to how much material can be captured from the accretion disc around a black hole, based on the power output, or luminosity of the infalling material. This limit is expressed in equation form by: L_edd =(4\u03C0GMmc)/\u03C3.\n"
              "20- A rearranged version of Newton's formula for universal gravitation:Here, `a' will be the difference in acceleration between two points separated by height `h', above a body of mass M and radius r.(a=2GMh/r^3)\n"
              "21- Falling Time. (T_fall = 15*10^-6(M/M_sun))\n"
              "22- Rotation of BH.The rotation of black holes is usually characterized by a number between zero and one, called a, which is calculated by the equation J is the angular momentum of the black hole, M is its mass, and G is Newton's gravitational Constant.(a=Jc/M^2G)\n"
              "23- Maximum Spin Rate Of BH.(F_max=16000*(M_sun/M))\n"
              "24- De Broglie Formula.(\u03BB = h/mv)\n"
              "25- Entropy Formula.(S=k.log(W)\n"
              "26- Entropy of BH.(S_BH = kA/4*l_p^2)\n"
              "27- Temperature Of BH.(T=k/2\u03C0)\n"
              "28- Temperature Of a Schwarzschild BH.(T= h_bar*c^3/8\u03C0GMk"
              "29- Wien's Law.(T=W/\u03BB_peak)\n"
              "30- Radius Of BH's Shadow.(R_shadow = sqrt(27)*GM/c^2 or 2.6*2GM/c^2)\n"
              "")
        formula_num = int(input("Please write your choice here (number of line): "))

        if formula_num == 0:
            key = False
        elif formula_num == 1:
            form_1()
        elif formula_num == 2:
            form_2()
        elif formula_num == 3:
            form_3()
        elif formula_num == 4:
            form_4()
        elif formula_num == 5:
            form_5()
        elif formula_num == 6:
            form_6()
        elif formula_num == 7:
            form_7()
        elif formula_num == 8:
            form_8()
        elif formula_num == 9:
            form_9()
        elif formula_num == 10:
            form_10()
        elif formula_num == 11:
            form_11()
        elif formula_num == 12:
            form_12()
        elif formula_num == 13:
            form_13()
        elif formula_num == 14:
            form_14()
        elif formula_num == 15:
            form_15()
        elif formula_num == 16:
            form_16()
        elif formula_num == 17:
            form_17()
        elif formula_num == 18:
            form_18()
        elif formula_num == 19:
            form_19()
        elif formula_num == 20:
            form_20()
        elif formula_num == 21:
            form_21()
        elif formula_num == 22:
            form_22()
        elif formula_num == 23:
            form_23()
        elif formula_num == 24:
            form_24()
        elif formula_num == 25:
            form_25()
        elif formula_num == 26:
            form_26()
        elif formula_num == 27:
            form_27()
        elif formula_num == 28:
            form_28()
        elif formula_num == 29:
            form_29()
        elif formula_num == 30:
            form_30()
        else:
            print("The entry value is invalid!\n")




main()