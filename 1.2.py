from sdtoolbox.postshock import CJspeed,PostShock_fr,PostShock_eq
from sdtoolbox.thermo import soundspeed_eq
from sdtoolbox.reflections import reflected_eq
import PySimpleGUI as sg
import cantera as ct
import numpy as np
import datetime
import sys
import os
a1=['aramco2.cti', 'Blanquart2018.cti', 'Burke2012.cti', 'Davis2005.cti', 'ffcm1.cti', 'gri30.cti', 'gri30_highT.cti', 'gri30_ion.cti', 'hexaneFull.cti', 'hexanePartial.cti', 'hexaneReduced.cti', 'Keromnes2013.cti', 'Li2015.cti', 'Mevel2015.cti', 'Mevel2018.cti', 'nDodecane_Reitz.cti', 'ptcombust.cti']
a2=['h2br2.cti', 'h2o2.cti', 'Hong2011.cti', 'Mevel2017.cti', 'sandiego20161214_H2only.cti']
b1=[['O2', 'N2', 'CH4/C2H6/C3H8/C2H4/C2H2/C2H5OH/C3H6/CH2O/CH3CHO/CH3COCH3/CH3OH'], ['AR', 'O2', 'N2', 'H2/CH4/C2H6/C3H8/C2H4/C3H4/C2H2/C2H5OH/C3H6/CH2O/CH3CHO/CH3COCH3/CH3OH/C3H8/C4H2/C4H6/C5H6/C6H6/C7H8/C8H8/C8H10/C11H10/C6H14/C7H16/C8H18/C12H26'], ['H2', 'H2O', 'O2', 'N2', 'CO', 'CO2', 'AR'], ['H2', 'AR', 'N2', 'CO'], ['AR', 'HE', 'N2', 'H2', 'CO', 'C'], ['H2', 'O2', 'CH4', 'N2', 'AR'], ['H2', 'O2', 'CH4', 'N2', 'AR'], ['H2', 'O2', 'CH4', 'N2', 'AR'], ['C6H14', 'O2', 'N2', 'AR'], ['N2','O2', 'C6H14', 'AR'], ['C6H14', 'O2', 'N2', 'AR'], ['H2', 'CO', 'O2', 'N2', 'AR', 'HE'], ['AR', 'N2', 'H2', 'O2', 'CO', 'HE'], ['H2', 'O2', 'NO2', 'N2', 'CH4/C2H6/C2H4/C2H2/C3H8/C3H6/C3H4'], ['CH3CHO', 'O2', 'N2', 'AR'], ['C12H26', 'O2', 'N2', 'CO2'], ['CH4', 'O2', 'N2', 'AR']]
b2=[['H2', 'BR2'], ['H2', 'O2', 'AR'], ['H2', 'O2', 'N2'], ['H2', 'O2', 'N2'], ['N2', 'H2', 'O2']]
a3=['aramco2.cti','Mevel2017.cti','Mevel2015.cti']
b3=[['CH4/C2H6/C3H8/C2H4/C2H2/C2H5OH/C3H6/CH2O/CH3CHO/CH3COCH3/CH3OH'],['H2'],[ 'H2','NO2', 'CH4/C2H6/C2H4/C2H2/C3H8/C3H6/C3H4']]
#dictm字典中储存每个反应物的质量
dictm={'H2':2,'CH3CHO':44,'CH4':16,'C2H6':30,'C3H8':44,'C2H4':28,'C2H2':26,'C2H50H':46,'C3H6':42,'CH20':30,'CH3CH0':44,'CH3C0CH3':58,'CH3OH':32,'C3H4':40}
#dictn字典储存每个反应物完全反应所需的氧气的物质的量
dictn={'H2':0.5,'CH3CHO':2.5,'CH4':2,'C2H6':3.5,'C3H8':5,'C2H4':3,'C2H2':2.5,'C2H50H':3,'C3H6':4.5,'CH20':1,'CH3CH0':2.5,'CH3C0CH3':4,'CH3OH':1.5,'C3H4':4}

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
def atmopara_cal(h):
    from math import  e
    if h<=11000:
        t0=(288.15-0.0065*h)#温度K
        p0=1.01325*10**5*(1-0.225577*10**(-4)*h)**5.25588#压强PA
        density0=1.225*(1-0.225577*10**(-4)*h)**4.25588#密度KG/M^3
    elif h>11000 and h<=20000:
        t0=216.5
        p0=2.263204*10**(-4)*e**(-1.576885*10**(-4)*(h-11000))
        density0=0.3639176*e**(-1.576885*10**(-4)*(h-11000))
    elif h>20000 and h<=32000:
        t0=216.5+0.001*(h-20000)
        p0=5.474849*10**3*(1+4.61574*10^(-6)*(h-20000))**(-34.16322)
        density0=8.803471*10**(-2)*(1+4.61574*10**(-6)*(h-20000))**(-34.16322)
    else:
        print("error in")
    return t0,p0,density0

#用于计算比冲，其中dp，A,L需要输入，Ucj为反应物燃烧的CJ速度，可由前部分计算给出

def analysis1(dp,A,L,Ucj):
    v=A*L
    i=4.3*dp*v/Ucj
    return i

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
def pda2 (yqsumn):
    for i in range(lens):
         try:
             a=float(valuesb['fywn{}'.format(i)])
             b=float(valuesb['c0'])
             d=float(valuesb['T1'])
             c=float(valuesb['P1'])
         except:
             sg.Popup('请输入正确的反应组分和初始条件！')
             return 0
    global q
    q = ''
    if pd2 == 1:
        b3[p].pop(-1)
        b3[p].append(valuesb['fywmulty'])
    for i in range(lens):
        q += b3[p][i] + ':' + valuesb['fywn{}'.format(i)] + ' '
    q+='O2'+':'+ '{}'.format(yqsumn) + ' '
    q+='N2'+':'+ '{}'.format(yqsumn*0.78/0.21) + ' '
    global T1
    T1 = float(valuesb['T1'])
    global c0
    c0 = float(valuesb['c0'])
    global P1
    P1 = float(valuesb['P1'])
    global mech
    mech = str
    global gas1
    gas1 = ct.Solution(mech)
    gas1.TPX = T1, P1, q
    global rho1
    rho1 = gas1.density
    print(q)
def pda ():
    for i in range(lens):
         try:
             a=float(values3['fywn{}'.format(i)])
             b=float(values3['T1'])
             c=float(values3['P1'])
         except:
             sg.Popup('请输入正确的反应组分和初始条件！')
             return 0
    if pd4!=1:
        global q
        q = ''
    sg.Popup('计算过程中请不要点击鼠标，否则程序可能未响应！\n点击确定开始计算')
    if pd == 1:
        if pd2 == 1:
            b1[p].pop(-1)
            b1[p].append(values3['fywmulty'])
        if pd4!=1:
            for i in range(lens):
                q += b1[p][i] + ':' + values3['fywn{}'.format(i)] + ' '
    elif pd == 2:
        if pd2 == 1:
            b2[p].pop(-1)
            b2[p].append(values3['fywmulty'])
        if pd4!= 1:
            for i in range(lens):
                q += b2[p][i] + ':' + values3['fywn{}'.format(i)] + ' '
    global T1
    T1 = float(values3['T1'])
    global P1
    P1 = float(values3['P1'])
    global mech
    mech = values2['fytype']
    global gas1
    gas1 = ct.Solution(mech)
    gas1.TPX = T1, P1, q
    global rho1
    rho1 = gas1.density
#   print(q)
# set initial state, composition, and gas object
q=''
h=''
A=''
L=''
Ma=''
f=''
c0=''
pd2=0
pd4=0
pd6=0
sg.change_look_and_feel('Dark blue 3')
lay0=[[sg.Text('    热力学计算程序',font=("times new roman", 16),size=(20,2))],
      [sg.Text('请选择计算类别',font=("times new roman", 16),size=(20,1))],
      [sg.Button('脉冲爆震发动机性能计算',font=("times new roman", 12),size=(30,2))],
      [sg.Button('爆震热力参数计算',font=("times new roman", 12),size=(30,2))],
    ]
win0=sg.Window(title='爆震热力参数计算',size=(300,250), default_element_size=(40, 1), grab_anywhere=True).Layout(lay0)
win1_active=False
wina_active=False
while 1:
    event0, values0 = win0.read()
    if event0 in (None, '退出'):
        sys.exit()
    if event0=='脉冲爆震发动机性能计算':
        wina_active=True
        win0.Hide()
        laya=[[sg.Text('请选择反应物成分：',font=("times new roman", 16),size=(20,1))],
            [sg.Radio('1.H2 O2 N2','RADIO1', default=True),sg.Radio('2.CH O2 N2','RADIO1', default=True),sg.Radio('3.CH H2 O2 N2','RADIO1', default=True)],
            [sg.Button('确定'),sg.Button('返回')]
            ]
        wina=sg.Window(title='脉冲爆震发动机性能计算',size=(400,120), default_element_size=(40, 1), grab_anywhere=True).Layout(laya)
        winb_active=False
        winc_active=False
        while 1:
           eventa, valuesa = wina.read()
           if eventa in (None, '退出'):
                sys.exit()
           if eventa == '返回':
               wina.close()
               wina_active = False
               win0.UnHide()
               break
           if eventa=='确定':
               menu_def = [['文件', ['打开', '保存']], ['编辑', ['默认','清除']], ['帮助', '关于...'], ]
               layb = [[sg.Menu(menu_def, tearoff=True)]]
               if valuesa[0]==True:
                   str='Mevel2017.cti'
                   p=a3.index('Mevel2017.cti')
                   lens=len(b3[p])
                   layb+=[[sg.Text('默认反应机理为:', font=("times new roman", 15),size=(40, 1),)],
                         [sg.Text('Mevel2017.cti', font=("times new roman", 18), size=(40, 1),text_color=('red'))]]
                   for i in range(lens):
                      if b3[p][i].count('/')!=0:
                          code=b3[p][i].split('/')
                          layb+=[
                              [sg.InputOptionMenu(code,size=(15,1),key='fywmulty'),sg.Text('',size=(8,1)),sg.Input(size=(35,1),key=('fywn{}'.format(i)))]
                          ]
                          pd2=1
                          break
                      layb+=[[sg.Text(b3[p][i], font=("times new roman", 15),size=(20, 1)),sg.Input(size=(35,1),key=('fywn{}'.format(i)))]]
               if valuesa[1]==True:
                   str='aramco2.cti'
                   p=a3.index('aramco2.cti')
                   lens=len(b3[p])
                   layb+=[[sg.Text('默认反应机理为:', font=("times new roman", 15),size=(40, 1),)],
                         [sg.Text('aramco2.cti', font=("times new roman", 18), size=(40, 1),text_color=('red'))]]
                   for i in range(lens):
                      if b3[p][i].count('/')!=0:
                          code=b3[p][i].split('/')
                          layb+=[
                              [sg.InputOptionMenu(code,size=(15,1),key='fywmulty'),sg.Text('',size=(8,1)),sg.Input(size=(35,1),key=('fywn{}'.format(i)))]
                          ]
                          pd2=1
                          break
                      layb+=[[sg.Text(b3[p][i], font=("times new roman", 15),size=(20, 1)),sg.Input(size=(35,1),key=('fywn{}'.format(i)))]]
               if valuesa[2]==True:
                   str='Mevel2015.cti'
                   p=a3.index('Mevel2015.cti')
                   lens=len(b3[p])
                   layb+=[[sg.Text('默认反应机理为:', font=("times new roman", 15),size=(40, 1),)],
                         [sg.Text('Mevel2015.cti', font=("times new roman", 18), size=(40, 1),text_color=('red'))]]
                   for i in range(lens):
                      if b3[p][i].count('/')!=0:
                          code=b3[p][i].split('/')
                          layb+=[
                              [sg.InputOptionMenu(code,size=(15,1),key='fywmulty'),sg.Text('',size=(8,1)),sg.Input(size=(35,1),key=('fywn{}'.format(i)))]
                          ]
                          pd2=1
                          break
                      layb+=[[sg.Text(b3[p][i], font=("times new roman", 15),size=(20, 1)),sg.Input(size=(35,1),key=('fywn{}'.format(i)))]]
               layb += [
               [sg.Text('海拔高度H/m', font=("times new roman", 15), size=(20, 1)),sg.Input(size=(35, 1), key=('h'))],
               [sg.Text('爆震室横截面积A/(m*m)', font=("times new roman", 15), size=(20, 1)),sg.Input(size=(35, 1), key=('A'))],
               [sg.Text('爆震室长度L/m', font=("times new roman", 15), size=(20, 1)), sg.Input(size=(35, 1),key=('L'))],
               [sg.Text('飞行马赫数Ma', font=("times new roman", 15), size=(20, 1)),sg.Input(size=(35, 1), key=('Ma'))],
               [sg.Text('爆震频率f', font=("times new roman", 15), size=(20, 1)),sg.Input(size=(35, 1), key=('f'))],
               [sg.Text('当量比', font=("times new roman", 15),size=(20, 1)),sg.Input(size=(35,1),key=('c0'))],
               [sg.Text('初始温度(K)', font=("times new roman", 15), size=(20, 1)), sg.Input(key='T1',size=(35,1))],
               [sg.Text('初始压强(Pa)', font=("times new roman", 15), size=(20, 1)), sg.Input(key='P1',size=(35,1))],
               [sg.Multiline('此处预览根据当量比计算出的参数量...', font=("times new roman", 15), size=(30, 5), key='out1'),sg.Button('预览')],
               [sg.Multiline('结果输出', font=("times new roman", 15), size=(80, 7), key='out2')],
               [sg.Button('确定'), sg.Button('返回')]
               ]
               winb_active=True
               wina.Hide()       
               winb = sg.Window(title='爆震热力参数计算',size=(500,800), default_element_size=(40, 10), grab_anywhere=True).Layout(layb)
               while 1:
                   eventb, valuesb = winb.read()
                   if eventb in (None, '退出'):
                       sys.exit()
                   if eventb == '默认':
                       winb['T1'].update(300)
                       winb['P1'].update(101300)
                       winb['h'].update(100)
                       winb['A'].update(100)
                       winb['L'].update(100)
                       winb['Ma'].update(1)
                       winb['f'].update(1)
                       winb['c0'].update(1)
                       for i in range(lens):
                           winb['fywn{}'.format(i)].update(1)
                   if eventb == '清除':
                       winb['h'].update('')
                       winb['A'].update('')
                       winb['L'].update('')
                       winb['Ma'].update('')
                       winb['f'].update('')
                       winb['c0'].update('')
                       winb['c0'].update('')
                       winb['T1'].update('')
                       winb['P1'].update('')
                       for i in range(lens):
                           winb['fywn{}'.format(i)].update('')
                   if eventb == '关于...':
                       sg.Popup('About this program', 'Version 3.0', 'PySimpleGUI rocks...','By:Solo dance')
                   if eventb == '打开':
                       filename = sg.PopupGetFile('file to open', no_window=True)
                       os.system(filename)
                   if eventb == '保存':
                       sg.Popup('选择一个txt文档，数据将会保存在里面')
                       savename = sg.PopupGetFile('file to open', no_window=True)
                       if savename!='':
                           with open(savename, 'a') as f:
                                 f.write(valuesb['out2'])
                   if eventb == '返回':
                       winb.close()
                       winb_active = False
                       wina.UnHide()
                       break
                   if eventb == '预览':
                        pda2(1)
                        fywsumn=0 #反应完全所需氧气的物质的量
                        fywsumm=0 #反应物总质量
                        q1=q.split()
                        for i in range(len(q1)):
                            q2=q1[i].split(':')
                            if q2[0] in dictm:
                                fywsumn+=dictn[q2[0]]*float(q2[1])
                                fywsumm+=dictm[q2[0]]*float(q2[1])
                        yqsumn=fywsumn/c0
                        yqsumm=yqsumn*5*28
                        c1=fywsumm/yqsumm#质量比
                        pda2(yqsumn)
                        winb['out1'].update('02:{}\nN2:{}'.format(yqsumn,yqsumn*0.78/0.21))
                   if eventb == '确定':
                        try:
                            pda2(1)
                            fywsumn=0 #反应完全所需氧气的物质的量
                            fywsumm=0 #反应物总质量
                            q1=q.split()
                            for i in range(len(q1)):
                                q2=q1[i].split(':')
                                if q2[0] in dictm:
                                    fywsumn+=dictn[q2[0]]*float(q2[1])
                                    fywsumm+=dictm[q2[0]]*float(q2[1])
                            yqsumn=fywsumn/c0
                            yqsumm=yqsumn*5*28
                            c1=fywsumm/yqsumm#质量比
                            pda2(yqsumn)


                            for i in ['h','A','L','Ma','f','c0']:
                                try:
                                    valuesb[i]=float(valuesb[i])
                                except:
                                    sg.Popup('请输入正确的反应组分和初始条件！')
                                    break
                            h = float(valuesb['h'])
                            A = float(valuesb['A'])
                            L = float(valuesb['L'])
                            Ma = float(valuesb['Ma'])
                            f = float(valuesb['f'])
                            c0 = float(valuesb['c0'])
                            lay4 = [[sg.Text('计算较慢，请耐心等待！')],
                                    [sg.ProgressBar(100, orientation='h', size=(20, 20), key=('progressbar'))]
                                    ]
                            win4 = sg.Window('计算中...', default_element_size=(40, 10), grab_anywhere=True,
                                            disable_close=True).Layout(
                                lay4)
                            progress_bar = win4['progressbar']
                            event4, values4 = win4.read(timeout=10)
                            for i in range(0, 5):
                                event4, values4 = win4.read(timeout=100)
                                progress_bar.UpdateBar(i + 1)
                            for i in range(5, 10):
                                event4, values4 = win4.read(timeout=100)
                                progress_bar.UpdateBar(i + 1)
                            # Find CJ speed
                            cj_speed = CJspeed(P1, T1, q, mech)
                            #  print("反应的CJ速度是",cj_speed,"m/s")
                            # calculate the CJ speed and find the VN point
                            # find the VN point
                            gas = PostShock_fr(cj_speed, P1, T1, q, mech)
                            gas_cjstate = PostShock_eq(cj_speed, P1, T1, q, mech)
                            denref = gas.density  # 密度
                            Tref = gas.T  # VN温度
                            Pref = gas.P  # VN压强
                            uvn = cj_speed * rho1 / denref
                            gamma = gas.cp / gas.cv  # 比热比
                            for i in range(10, 30):
                                event4, values4 = win4.read(timeout=100)
                                progress_bar.UpdateBar(i + 1)
                            #  print("Tvn", Tref)
                            #  print ("Dvn", denref)
                            #  print ("gamma", gamma)
                            # get the nominal desired  ignition delay and reaction time
                            [tignition, texplosion] = characteristictimes(gas)
                            #  print ("tignition", tignition )
                            #  print ("texplosion", texplosion )

                            # 活化能
                            # perturb the temperature by +-perturbation to obtain the activation temperature Ta (in K)
                            # in ti = A rho^(-n) Exp (Ta/T)
                            # from the change in ignition delay, while keeping the density constant
                            perturbation = 0.1
                            Tplus = Tref * (1.0 + perturbation)
                            Tminus = Tref * (1.0 - perturbation)
                            gas.TDX = Tplus, denref, q
                            [tignitionplus, texplosionplus] = characteristictimes(gas)
                            gas.TDX = Tminus, denref, q
                            [tignitionminus, texplosionminus] = characteristictimes(gas)
                            activationtemperature = -(Tref / tignition) * (tignitionplus - tignitionminus) / (
                                    2.0 * perturbation)
                            for i in range(30, 70):
                                event4, values4 = win4.read(timeout=20)
                                progress_bar.UpdateBar(i + 1)
                            #  print ("activation temperature", activationtemperature)
                            # M_CJ  Ma=(r*R/M*T)**0.5  R——摩尔气体常数（R=8.314J/mol•K）；
                            #  r——比热容之比（气体定压比热容与定容比热容之比）；
                            #  M——平均摩尔质量；
                            #  T——气体的开氏温度。
                            R = ct.gas_constant
                            a = (gamma * R * T1 / (16 + 1 + 12)) ** 0.5
                            Mcj = cj_speed / a
                            Pcj=gas_cjstate.P
                            #  print ("Mcj" , Mcj)
                            # heat releas  Q/(R*T0) 放出热量
                            Q = (Mcj - 1 / Mcj) * (Mcj - 1 / Mcj) * gamma / 2 / (gamma * gamma - 1)
                            for i in range(70, 100):
                                event4, values4 = win4.read(timeout=10)
                                progress_bar.UpdateBar(i + 1)
                            #  print("Q/RT0", Q)
                            # cov = gas1.viscosity
                            #  print("粘性系数为",cov)
                            # d=gas.density
                            # cov1=gas1.viscosity
                            # d1=gas1.density
                            #  print("导热系数为",heatd)




                            # c1为燃料和空气的质量比例，根据选择反应机理和初始输入组分的不同给出
                            t0, p0, density0 = atmopara_cal(h)  # h为发动机工作时的飞行高度，需要输入给出
                            c0 = c0 = (1.4 * p0 / density0) ** 0.5  # c0为当地声速，1.4为比热比
                            u0 = Ma * c0  # Ma为飞行马赫数
                            ma = A * u0 * density0  # ma为空气流量
                            mfill = ma * c1  # c1为燃料和空气的质量比例，根据选择反应机理和初始输入组分的不同给出,mfill为燃料填充质量流量
                            c = 346  # 室温下声速
                            ##以下对p3进行计算
                            c3 = cj_speed / 2
                            c2 = 2 / (gamma + 1) * (c3 + (gamma - 1) / 2 * (cj_speed + c))
                            p3 = Pcj * (c3 / c2) ** (2 * gamma / (gamma - 1))  # p2为cj压力，是上部分输出的Pcj
                            dp = p3 - p0

                            i = analysis1(dp, A, L, cj_speed)  # 其中dp=p3-p0，A,L需要输入，Ucj为反应物燃烧的CJ速度，可由前部分计算给出

                            Frn = 0  # 壁面摩擦系数

                            if Ma < 1:

                                sita_i = 1
                            else:
                                k = 1.4
                                sita_i = ((k + 1) * Ma ** 2 / 2 / (1 + (k - 1) / 2 * Ma ** 2)) ** (k / (k - 1)) * (
                                            2 * k * Ma ** 2 / (k + 1) - (k - 1) / (k + 1)) ** (-1 / (k - 1))

                            sita_c = 1 - 4.468 * Frn * mfill ** 2.16  # 摩擦总压损失
                            P0 = p0 + 0.5 * density0 * u0 ** 2  # 进气端总压
                            P_ = P0 * sita_c * sita_i  # 出口总压
                            pe = 2 * P_ - (u0 ** 2 * density0 + p0 + c1 * u0 * density0 * c)  # 出口压强
                            ue = 2 * (u0 ** 2 * density0 + p0 + c1 * u0 * density0 * c - P_) / (
                                        u0 * density0 + u0 * density0 * c1)  # 出口速度
                            F = i * f + ma * (ue - u0) + (pe - p0) * A  # u0为发动机工作时的飞行速度，需要输入给出；f为爆震频率，需要输入给出。
                            winb['out2'].update('计算结果为：\n  直管爆震室单个工作周期冲量:{}\n 平均标准推力:{}\n  '.format(i, F))
                            win4.close()
                        except:
                            sg.Popup('反应出错！请检查反应物和反应机理！')





    if event0=='爆震热力参数计算':    #此处与1.1版本相同
        win1_active=True
        win0.Hide()
        lay1=[[sg.Text('爆震热力参数计算',font=("times new roman", 16),size=(20,2))],
              [sg.Text('请选择反应模型',font=("times new roman", 16),size=(20,1))],
              [sg.Button('碳氢空气反应模型',size=(15,1)),sg.Button('氢空气反应模型',size=(15,1))],
              [sg.Button('返回')]
              ]
        win1=sg.Window(title='爆震热力参数计算',size=(300,150), default_element_size=(40, 1), grab_anywhere=True).Layout(lay1)
        win2_active=False
        win3_active=False
        win5_active=False
        win6_active=False
        while 1:
           event1, values1 = win1.read()
           if event1 in (None, '退出'):
                sys.exit()
           if event1 == '返回':
               win1.close()
               win1_active = False
               win0.UnHide()
               break
           if event1=='氢空气反应模型':
               pd=2
               temp = pd
               win2_active=True
               win1.Hide()
               lay2 = [[sg.Text('爆震热力参数计算', font=("times new roman", 16), size=(20, 2))],
                       [sg.Text('您选择的反应模型为：', font=("times new roman", 16), size=(20, 1))],
                       [sg.Text('氢空气反应模型', font=("times new roman", 12), size=(20, 1),text_color=('red'),key='fymod')],
                       [sg.Text('请在下面选择反应机理', font=("times new roman", 16), size=(20, 1))],
                       [sg.InputOptionMenu(('h2br2.cti','h2o2.cti','Hong2011.cti','Mevel2017.cti','sandiego20161214_H2only.cti'),key='fytype')],
                       [sg.Button('确定'),sg.Button('返回')]
                       ]
               str='氢空气反应模型'
               win2 = sg.Window(title='爆震热力参数计算',size=(300,250), default_element_size=(40, 10), grab_anywhere=True).Layout(lay2)
           if event1=='碳氢空气反应模型':
               pd=1
               temp = pd
               win2_active=True
               win1.Hide()
               lay2 = [[sg.Text('爆震热力参数计算', font=("times new roman", 16), size=(20, 2))],
                       [sg.Text('您选择的反应模型为：', font=("times new roman", 16), size=(20, 1))],
                       [sg.Text('碳氢空气反应模型', font=("times new roman", 12), size=(20, 1),text_color=('red'))],
                       [sg.Text('请在下面选择反应机理', font=("times new roman", 16), size=(20, 1))],
                       [sg.InputOptionMenu(('aramco2.cti','Blanquart2018.cti','Burke2012.cti','Davis2005.cti','ffcm1.cti','gri30.cti','gri30_highT.cti','gri30_ion.cti','hexaneFull.cti','hexanePartial.cti','hexaneReduced.cti','Keromnes2013.cti','Li2015.cti','methane_pox_on_pt.cti','Mevel2015.cti','Mevel2018.cti','nDodecane_Reitz.cti','ptcombust.cti','sandiego20161214.cti'),key='fytype')],
                       [sg.Button('确定'),sg.Button('返回')]
                       ]
               str='碳氢空气反应模型'
               win2 = sg.Window(title='爆震热力参数计算',size=(300,270), default_element_size=(40, 10), grab_anywhere=True).Layout(lay2)
           while 1:
               event2, values2 = win2.read()
               if event2 in (None, '退出'):
                   sys.exit()
               if event2 == '返回':
                   win2.close()
                   win2_active = False
                   win1.UnHide()
                   break
               if event2 == '确定':
                   win3_active = True
                   win2.Hide()
                   menu_def = [['文件', ['打开', '保存']], ['编辑', ['默认','清除']], ['帮助', '关于...'], ]
                   lay3 = [[sg.Menu(menu_def, tearoff=True)],
                           [sg.Text('选择的反应机理为:', font=("times new roman", 15),size=(60, 1),)],
                           [sg.Text(values2['fytype'], font=("times new roman", 18), size=(40, 1),text_color=('red'))]
                   ]
                   if str=='碳氢空气反应模型':
                       p=a1.index(values2['fytype'])
                       lens=len(b1[p])
                       for i in range(lens):
                           pd2=0
                           if b1[p][i].count('/')!=0:
                               pd2=1
                               code=b1[p][i].split('/')
                               lay3+=[
                                 [sg.InputOptionMenu(code,size=(15,1),key='fywmulty'),sg.Text('',size=(8, 1)),sg.Input(size=(35,1),key=('fywn{}'.format(i)))]
                               ]
                               break
                           lay3+=[
                               [sg.Text(b1[p][i], font=("times new roman", 15),size=(20, 1)),sg.Input(size=(35,1),key=('fywn{}'.format(i)))]
                           ]
                   elif str == '氢空气反应模型':
                       p=a2.index(values2['fytype'])
                       lens = len(b2[p])
                       for i in range(lens):
                           pd2=0
                           if b1[p][i].count('/')!=0:
                               pd2=1
                               code=b1[p][i].split('/')
                               lay3+=[
                                   [sg.InputOptionMenu(code,size=(15,1),key='fywmulty'),sg.Input(size=(35,1),key=('fywn{}'.format(i)))]
                               ]
                               break
                           lay3+=[
                               [sg.Text(b2[p][i], font=("times new roman", 15),size=(20, 1)),sg.Input(size=(35,1),key=('fywn{}'.format(i)))]
                           ]
                   lay3 += [
                       [sg.Text('初始温度(K)', font=("times new roman", 15), size=(20, 1)), sg.Input(key='T1',size=(35,1))],
                       [sg.Text('初始压强(Pa)', font=("times new roman", 15), size=(20, 1)), sg.Input(key='P1',size=(35,1)),sg.Button('自定义...',size=(15,1))],
                       [sg.Multiline('结果输出', font=("times new roman", 15), size=(27, 8), key='out1'),sg.Multiline('其他系数输出', font=("times new roman", 15), size=(27,8), key='out2')],
                       [sg.Button('粘度系数'), sg.Button('质量扩散系数'), sg.Button('导热系数')],
                       [sg.Button('确定'), sg.Button('返回')],
                   ]
                   win3 = sg.Window('爆震热力参数计算', default_element_size=(40, 10), grab_anywhere=True).Layout(lay3)
               while 1:
                   event3,values3=win3.read()
                   if event3 in (None, '退出'):
                       sys.exit()
                   elif event3 == '自定义...':
                       win5_active=True
                       lay5=[
                           [sg.Text('已设定的q值为：'),sg.Text(q,size=(50,1),key='input')],
                           [sg.Text('q=', font=("times new roman", 15), size=(20, 1)), sg.Input(key='q',size=(35,1))],
                           [sg.Button('确定'),sg.Button('还原'),sg.Button('返回'),sg.Text('',size=(20,1)),sg.Button('帮助',size=(14,1))]
                       ]
                       win5 = sg.Window('自定义...', default_element_size=(40, 10), grab_anywhere=True).Layout(lay5)
                       while 1:
                           event5,values5=win5.read()
                           if event5 in (None, '返回'):
                               win5.close()
                               win5_active=False
                               if pd4==1:
                                   for i in range(lens):
                                       win3['fywn{}'.format(i)].update('0000')
                               break
                           elif event5=='还原':
                               pd4=0
                               q=''
                               sg.Popup('已还原q的数值',font=("times new roman", 15))
                               win5['input'].update(q)
                           elif event5=='确定':
                               pd3 = 1
                               q=values5['q']
                               m = q.split(' ')
                               mm = []
                               for i in range(len(m)):
                                   mm.append(m[i].split(':'))
                               for i in range(len(mm)):
                                   if len(mm[i]) != 2:
                                       pd3 = 0
                                       break
                                   if type(mm[i][-1]) != float:
                                       pd3 = 0
                                       break
                               if pd3==0:
                                   sg.Popup('输入q的格式不正确！\n请点击帮助获得提示')
                               elif pd3==1:
                                   pd4=1
                                   sg.Popup('已设定q的数值', font=("times new roman", 15))
                                   win5['input'].update(q)
                           elif event5=='帮助':
                               sg.Popup('输入的q为所有反应物与其组分的组合\n格式示例 \'q=C2H4:1 H2:1 O2:1\'\n输入结束后点击确定，原位置输入的组分全部置为0000(作废)，请不要随意更改！\n点击还原可还原自定义输入',title='帮助',font=("times new roman", 15))
                   elif event3 == '返回':
                       win3.close()
                       win3_active = False
                       win2.UnHide()
                       break
                   elif event3 == '粘度系数':
                       try:
                           pda()
                           gas1.transport_model = 'Multi'
                           cov = gas1.viscosity
                           win3['out2'].update('粘度系数为：{}'.format(cov))
                       except:
                           sg.Popup('反应出错！请检查反应物和反应机理！')
                   elif event3 == '质量扩散系数':
                       try:
                           pda()
                           # 质量扩散系数  扩散系数是指当浓度梯度为一个单位时，单位时间内通过单位面积的气体量，
                           massdiff = gas1.mix_diff_coeffs_mass
                           win3['out2'].update('质量扩散系数为：{}'.format(massdiff))
                       except:
                           sg.Popup('反应出错！请检查反应物和反应机理！')
                   elif event3 == '导热系数':
                       try:
                           pda()
                           # 导热系数  导热系数是指在稳定传热条件下，1m厚的材料，两侧表面的温差为1度（K，°C）,在1秒内，通过1平方米面积传递的热量，用λ表示，
                           # 单位为瓦/米·度，W/m·k（W/m·K,此处的K可用℃代替）。
                           heatd = gas1.thermal_conductivity
                           win3['out2'].update('导热系数为：{}'.format(heatd))
                       except:
                           sg.Popup('反应出错！请检查反应物和反应机理！')
                   elif event3 == '默认':
                       win3['T1'].update(300)
                       win3['P1'].update(101300)
                       for i in range(lens):
                           win3['fywn{}'.format(i)].update(1)
                   elif event3 == '清除':
                       win3['T1'].update('')
                       win3['P1'].update('')
                       for i in range(lens):
                           win3['fywn{}'.format(i)].update('')
                   elif event3 == '关于...':
                       sg.Popup('About this program', 'Version 3.0', 'PySimpleGUI rocks...','By:Solo dance')
                   elif event3 == '打开':
                       filename = sg.PopupGetFile('file to open', no_window=True)
                       os.system(filename)
                   elif event3 == '保存':
                       sg.Popup('选择一个txt文档，数据将会保存在里面')
                       savename = sg.PopupGetFile('file to open', no_window=True)
                       if savename!='':
                           with open(savename, 'a') as f:
                                 f.write(values3['out1'])
                                 f.write(values3['out2'])
                   elif event3 == '确定':
                       try:
                           pda()
                           lay4 = [[sg.Text('计算较慢，请耐心等待！')],
                                   [sg.ProgressBar(100, orientation='h', size=(20, 20), key=('progressbar'))]
                                   ]
                           win4 = sg.Window('计算中...', default_element_size=(40, 10), grab_anywhere=True,
                                            disable_close=True).Layout(
                               lay4)
                           progress_bar = win4['progressbar']
                           event4, values4 = win4.read(timeout=10)
                           for i in range(0, 5):
                               event4, values4 = win4.read(timeout=100)
                               progress_bar.UpdateBar(i + 1)
                           for i in range(5, 10):
                               event4, values4 = win4.read(timeout=100)
                               progress_bar.UpdateBar(i + 1)
                           # Find CJ speed
                           cj_speed = CJspeed(P1, T1, q, mech)
                           #  print("反应的CJ速度是",cj_speed,"m/s")
                           # calculate the CJ speed and find the VN point
                           # find the VN point
                           gas = PostShock_fr(cj_speed, P1, T1, q, mech)
                           gas_cjstate = PostShock_eq(cj_speed, P1, T1, q, mech)
                           denref = gas.density  # 密度
                           Tref = gas.T  # VN温度
                           Pref = gas.P  # VN压强
                           uvn = cj_speed * rho1 / denref
                           gamma = gas.cp / gas.cv  # 比热比
                           for i in range(10, 30):
                               event4, values4 = win4.read(timeout=100)
                               progress_bar.UpdateBar(i + 1)
                           #  print("Tvn", Tref)
                           #  print ("Dvn", denref)
                           #  print ("gamma", gamma)
                           # get the nominal desired  ignition delay and reaction time
                           [tignition, texplosion] = characteristictimes(gas)
                           #  print ("tignition", tignition )
                           #  print ("texplosion", texplosion )

                           # 活化能
                           # perturb the temperature by +-perturbation to obtain the activation temperature Ta (in K)
                           # in ti = A rho^(-n) Exp (Ta/T)
                           # from the change in ignition delay, while keeping the density constant
                           perturbation = 0.1
                           Tplus = Tref * (1.0 + perturbation)
                           Tminus = Tref * (1.0 - perturbation)
                           gas.TDX = Tplus, denref, q
                           [tignitionplus, texplosionplus] = characteristictimes(gas)
                           gas.TDX = Tminus, denref, q
                           [tignitionminus, texplosionminus] = characteristictimes(gas)
                           activationtemperature = -(Tref / tignition) * (tignitionplus - tignitionminus) / (
                                   2.0 * perturbation)
                           for i in range(30, 70):
                               event4, values4 = win4.read(timeout=20)
                               progress_bar.UpdateBar(i + 1)
                           #  print ("activation temperature", activationtemperature)
                           # M_CJ  Ma=(r*R/M*T)**0.5  R——摩尔气体常数（R=8.314J/mol•K）；
                           #  r——比热容之比（气体定压比热容与定容比热容之比）；
                           #  M——平均摩尔质量；
                           #  T——气体的开氏温度。
                           R = ct.gas_constant
                           a = (gamma * R * T1 / (16 + 1 + 12)) ** 0.5
                           Mcj = cj_speed / a
                           Pcj=gas_cjstate.P
                           #  print ("Mcj" , Mcj)
                           # heat releas  Q/(R*T0) 放出热量
                           Q = (Mcj - 1 / Mcj) * (Mcj - 1 / Mcj) * gamma / 2 / (gamma * gamma - 1)
                           for i in range(70, 100):
                               event4, values4 = win4.read(timeout=10)
                               progress_bar.UpdateBar(i + 1)
                           #  print("Q/RT0", Q)
                           # cov = gas1.viscosity
                           #  print("粘性系数为",cov)
                           # d=gas.density
                           # cov1=gas1.viscosity
                           # d1=gas1.density
                           #  print("导热系数为",heatd)
                           win3['out1'].update(
                               '计算结果为：\n  CJ速度:{}\n CJ压力:{}\n Tvn:{}\n Dvn:{}\n  gamma:{}\n  tignition:{}\n  texplosion:{}\n  activation temperature:{}\n  Mcj:{}\n  Q/RTO:{}\n  '.format(
                                   cj_speed,Pcj,Tref, denref, gamma, tignition, texplosion, activationtemperature, Mcj, Q
                               ))
                           win4.close()
                       except:
                           sg.Popup('反应出错！请检查反应物和反应机理！')






