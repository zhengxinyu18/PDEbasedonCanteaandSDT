# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

from sdtoolbox.postshock import CJspeed,PostShock_fr
import PySimpleGUI as sg
import cantera as ct
import os

def getthermicity(gas):
    """
    Returns the thermicity = sum ( (w/wi-hsi/(cp*T))*dyidt )
    SYNTAX: thermicity = getthermicity(gas)
    """
    j = gas.n_species  #number of species  =   53
    thermicity=0.0
    for i in range(j):
        w=gas.mean_molecular_weight  #mean molecular weight of the mixture
        wi=gas.molecular_weights[i]  #molecular weight of the ith specie
        R = ct.gas_constant          #universal gas constant
        T =gas.T
        hsi=gas.standard_enthalpies_RT[i]*R*T/wi #specific enthalpy of the ith specie
        cp=gas.cp_mass               #mixture averaged cp, frozen
        rho=gas.density
        dyidt = gas.net_production_rates[i]*wi/rho  #dYi/dt: rate of change of the mass fraction of the ith specie
        thermicity =  thermicity + (w/wi-hsi/(cp*T))*dyidt
    return thermicity

def characteristictimes(gas):
    """
    Returns the ignition and reaction time from constant volume calculations
    [tignition, texplosion] = characteristictimes(gas)
    """
    combustor = ct.IdealGasReactor(gas)
    sim = ct.ReactorNet([combustor])
    tfinal = 2     #final time in the calculation
    tnow = 0.0
    told = 0.0
    thermicity = 0.0
    thermicitymax = 0.0
    tignition = 0.0
    texplosion = 0.0
    while tnow < tfinal:
        tnow = sim.step()
        thermicity=getthermicity(gas)
        told = tnow
        if thermicity > thermicitymax:
           thermicitymax = thermicity
           tignition = tnow
           texplosion = 1.0/thermicitymax
    del combustor
    del sim
    return [tignition, texplosion]

print(os.getcwd())
# set initial state, composition, and gas object
sg.change_look_and_feel('DarkAmber')
layout=[
    [sg.Text('请输入初始值:')],
    [sg.Text('t1',size=(15,1)),sg.Input(key='T1')],
    [sg.Text('p1',size=(15,1)),sg.Input(key='P1')],
    [sg.Text('ic2h2',size=(15,1)),sg.Input(key='iC2H2')],
    [sg.Text('io2',size=(15,1)),sg.Input(key='iO2')],
    [sg.Button('确定'),sg.Button('退出')]
]
window=sg.Window('技术创新',layout)
event, values = window.read()
while 1:
    if event in (None,'退出'):
        break
    elif event=='确定':
        sg.Popup('计算结果较慢，请耐心等待！')
        q='C2H2:'+values['iC2H2']+' O2:'+values['iO2']
        T1=int(values['T1'])
        P1=int(values['P1'])
        mech = '.\Mevel2015.cti'
        gas1 = ct.Solution(mech)

        gas1.TPX = T1,P1,q
        rho1 = gas1.density

# Find CJ speed
        cj_speed = CJspeed(P1, T1, q, mech)
#  print("反应的CJ速度是",cj_speed,"m/s")
# calculate the CJ speed and find the VN point

# find the VN point
        gas = PostShock_fr(cj_speed, P1, T1, q, mech)
        denref = gas.density #密度
        Tref = gas.T    #VN温度
        Pref = gas.P   #VN压强
        uvn = cj_speed * rho1 / denref
        gamma = gas.cp / gas.cv  #比热比

#  print("Tvn", Tref)
#  print ("Dvn", denref)
#  print ("gamma", gamma)
# get the nominal desired  ignition delay and reaction time
        [tignition, texplosion] = characteristictimes(gas)

#  print ("tignition", tignition )
#  print ("texplosion", texplosion )

#活化能
# perturb the temperature by +-perturbation to obtain the activation temperature Ta (in K)
# in ti = A rho^(-n) Exp (Ta/T)
# from the change in ignition delay, while keeping the density constant
        perturbation=0.1
        Tplus=Tref*(1.0+perturbation)
        Tminus=Tref*(1.0-perturbation)
        gas.TDX = Tplus, denref, q
        [tignitionplus, texplosionplus] = characteristictimes(gas)
        gas.TDX = Tminus, denref, q
        [tignitionminus, texplosionminus] = characteristictimes(gas)
        activationtemperature=-(Tref/tignition)*(tignitionplus-tignitionminus)/(2.0*perturbation)
#  print ("activation temperature", activationtemperature)

#M_CJ  Ma=(r*R/M*T)**0.5  R——摩尔气体常数（R=8.314J/mol•K）；
      #  r——比热容之比（气体定压比热容与定容比热容之比）；
      #  M——平均摩尔质量；
      #  T——气体的开氏温度。
        R = ct.gas_constant
        a = (gamma*R*T1/(16+1+12))**0.5
        Mcj = cj_speed/a
#  print ("Mcj" , Mcj)

# heat releas  Q/(R*T0) 放出热量
        Q=(Mcj-1/Mcj)*(Mcj-1/Mcj)*gamma/2/(gamma*gamma-1)
#  print("Q/RT0", Q)

#粘度系数  牛顿内摩擦定律的数学表达式为F=S*(du/dy)。而面积上的摩擦力(切应力)为F/S=μ(du/dy)。μ称为动力粘度或动力粘性系数。
# 把μ与流体密度P的比值称为运动粘度或运动粘性系数V。
        gas1.equilibrate
        cov=gas1.viscosity
#  print("粘性系数为",cov)
#d=gas.density
#cov1=gas1.viscosity
#d1=gas1.density
#  print("粘性系数为",cov/d)
#  print("粘性系数1为",cov1/d1)

#质量扩散系数  扩散系数是指当浓度梯度为一个单位时，单位时间内通过单位面积的气体量，
        massdiff=gas1.mix_diff_coeffs_mass
#print("质量扩散系数为",massdiff)

#导热系数  导热系数是指在稳定传热条件下，1m厚的材料，两侧表面的温差为1度（K，°C）,在1秒内，通过1平方米面积传递的热量，用λ表示，
# 单位为瓦/米·度，W/m·k（W/m·K,此处的K可用℃代替）。
        heatd=gas1.thermal_conductivity
#  print("导热系数为",heatd)
        sg.Popup('计算结果为：\n  CJ速度:{}\n  Tvn:{}\n  gamma:{}\n  tignition:{}\n  texplosion:{}\n  activation temperature:{}\n  Mcj:{}\n  Q/RTO:{}\n  质量扩散系数:{}\n  导热系数:{}\n  粘性系数：{}\n'.format(cj_speed,Tref,denref,gamma,tignition,texplosion,activationtemperature,Mcj,Q,massdiff,heatd,cov))
        event, values = window.read()
window.close()

