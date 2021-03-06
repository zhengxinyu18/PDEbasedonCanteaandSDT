
from sdtoolbox.postshock import CJspeed,PostShock_fr
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
def pda ():
    for i in range(lens):
        if values3['fywn{}'.format(i)] == '' or values3['T1'] == '' or values3['P1'] == '':
            sg.Popup('??????????????????????????????????????????')
            return 0
        if values3['fywn{}'.format(i)].isdigit() == False or values3['T1'].isdigit() == False or values3['P1'].isdigit() == False:
            sg.Popup('????????????????????????????????????????????????')
            return 0
    if pd4!=1:
        global q
        q = ''
    sg.Popup('?????????????????????????????????????????????????????????????????????\n????????????????????????')
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
    T1 = int(values3['T1'])
    global P1
    P1 = int(values3['P1'])
    global mech
    mech = values2['fytype']
    global gas1
    gas1 = ct.Solution(mech)
    gas1.TPX = T1, P1, q
    global rho1
    rho1 = gas1.density
# set initial state, composition, and gas object
q=''
pd4 = 0
sg.change_look_and_feel('Dark blue 3')
lay1=[[sg.Text('???????????????????????????\n????????????????????????',font=("times new roman", 16),size=(20,2))],
      [sg.Text('?????????????????????',font=("times new roman", 16),size=(20,1))],
      [sg.Button('????????????????????????',size=(15,1)),sg.Button('?????????????????????',size=(15,1))]
      ]
win1=sg.Window(title='???????????????????????????',size=(300,150), default_element_size=(40, 1), grab_anywhere=True).Layout(lay1)
win2_active=False
win3_active=False
win5_active=False
while 1:
   event1, values1 = win1.read()
   if event1 in (None, '??????'):
        sys.exit()
   if event1=='?????????????????????':
       pd=2
       temp = pd
       win2_active=True
       win1.Hide()
       lay2 = [[sg.Text('???????????????????????????\n????????????????????????', font=("times new roman", 16), size=(20, 2))],
               [sg.Text('??????????????????????????????', font=("times new roman", 16), size=(20, 1))],
               [sg.Text('?????????????????????', font=("times new roman", 12), size=(20, 1),text_color=('red'),key='fymod')],
               [sg.Text('??????????????????????????????', font=("times new roman", 16), size=(20, 1))],
               [sg.InputOptionMenu(('h2br2.cti','h2o2.cti','Hong2011.cti','Mevel2017.cti','sandiego20161214_H2only.cti'),key='fytype')],
               [sg.Button('??????'),sg.Button('??????')]
               ]
       str='?????????????????????'
       win2 = sg.Window(title='???????????????????????????',size=(300,250), default_element_size=(40, 10), grab_anywhere=True).Layout(lay2)
   if event1=='????????????????????????':
       pd=1
       temp = pd
       win2_active=True
       win1.Hide()
       lay2 = [[sg.Text('???????????????????????????\n????????????????????????', font=("times new roman", 16), size=(20, 2))],
               [sg.Text('??????????????????????????????', font=("times new roman", 16), size=(20, 1))],
               [sg.Text('????????????????????????', font=("times new roman", 12), size=(20, 1),text_color=('red'))],
               [sg.Text('??????????????????????????????', font=("times new roman", 16), size=(20, 1))],
               [sg.InputOptionMenu(('aramco2.cti','Blanquart2018.cti','Burke2012.cti','Davis2005.cti','ffcm1.cti','gri30.cti','gri30_highT.cti','gri30_ion.cti','hexaneFull.cti','hexanePartial.cti','hexaneReduced.cti','Keromnes2013.cti','Li2015.cti','methane_pox_on_pt.cti','Mevel2015.cti','Mevel2018.cti','nDodecane_Reitz.cti','ptcombust.cti','sandiego20161214.cti'),key='fytype')],
               [sg.Button('??????'),sg.Button('??????')]
               ]
       str='????????????????????????'
       win2 = sg.Window(title='???????????????????????????',size=(300,250), default_element_size=(40, 10), grab_anywhere=True).Layout(lay2)
   while 1:
       event2, values2 = win2.read()
       if event2 in (None, '??????'):
           sys.exit()
       if event2 == '??????':
           win2.close()
           win2_active = False
           win1.UnHide()
           break
       if event2 == '??????':
           win3_active = True
           win2.Hide()
           menu_def = [['??????', ['??????', '??????']], ['??????', ['??????','??????']], ['??????', '??????...'], ]
           lay3 = [[sg.Menu(menu_def, tearoff=True)],
                   [sg.Text('????????????????????????:', font=("times new roman", 15),size=(40, 1),)],
                   [sg.Text(values2['fytype'], font=("times new roman", 18), size=(40, 1),text_color=('red'))]
           ]
           if str=='????????????????????????':
               p=a1.index(values2['fytype'])
               lens=len(b1[p])
               for i in range(lens):
                   pd2=0
                   if b1[p][i].count('/')!=0:
                       pd2=1
                       code=b1[p][i].split('/')
                       lay3+=[
                           [sg.InputOptionMenu(code,size=(15,1),key='fywmulty'),sg.Text('',size=(8,1)),sg.Input(size=(35,1),key=('fywn{}'.format(i)))]
                       ]
                       break
                   lay3+=[
                       [sg.Text(b1[p][i], font=("times new roman", 15),size=(20, 1)),sg.Input(size=(35,1),key=('fywn{}'.format(i)))]
                   ]
           elif str == '?????????????????????':
               p=a2.index(values2['fytype'])
               lens = len(b2[p])
               for i in range(lens):
                   pd2=0
                   if b1[p][i].count('/')!=0:
                       pd2=1
                       code=b1[p][i].split('/')
                       lay3+=[
                           [sg.InputOptionMenu(code,size=(15,1),key='fywmulty'),sg.Text('',size=(8,1)),sg.Input(size=(35,1),key=('fywn{}'.format(i)))]
                       ]
                       break
                   lay3+=[
                       [sg.Text(b2[p][i], font=("times new roman", 15),size=(20, 1)),sg.Input(size=(35,1),key=('fywn{}'.format(i)))]
                   ]
           lay3 += [
               [sg.Text('????????????(K)', font=("times new roman", 15), size=(20, 1)), sg.Input(key='T1',size=(35,1))],
               [sg.Text('????????????(Pa)', font=("times new roman", 15), size=(20, 1)), sg.Input(key='P1',size=(35,1)),sg.Button('?????????...',size=(15,1))],
               [sg.Multiline('????????????', font=("times new roman", 15), size=(40, 8), key='out1'),sg.Multiline('??????????????????', font=("times new roman", 15), size=(40,8), key='out2')],
               [sg.Button('????????????'), sg.Button('??????????????????'), sg.Button('????????????')],
               [sg.Button('??????'), sg.Button('??????')],
           ]
           win3 = sg.Window('???????????????????????????', default_element_size=(40, 10), grab_anywhere=True).Layout(lay3)
       while 1:
           event3,values3=win3.read()
           if event3 in (None, '??????'):
               sys.exit()
           elif event3 == '?????????...':
               win5_active=True
               lay5=[
                   [sg.Text('????????????q?????????'),sg.Text(q,size=(50,1),key='input')],
                   [sg.Text('q=', font=("times new roman", 15), size=(20, 1)), sg.Input(key='q',size=(35,1))],
                   [sg.Button('??????'),sg.Button('??????'),sg.Button('??????'),sg.Text('',size=(20,1)),sg.Button('??????',size=(14,1))]
               ]
               win5 = sg.Window('?????????...', default_element_size=(40, 10), grab_anywhere=True).Layout(lay5)
               while 1:
                   event5,values5=win5.read()
                   if event5 in (None, '??????'):
                       win5.close()
                       win5_active=False
                       if pd4==1:
                           for i in range(lens):
                               win3['fywn{}'.format(i)].update('0000')
                       break
                   elif event5=='??????':
                       pd4=0
                       q=''
                       sg.Popup('?????????q?????????',font=("times new roman", 15))
                       win5['input'].update(q)
                   elif event5=='??????':
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
                           if mm[i][-1].isdigit() == False:
                               pd3 = 0
                               break
                       if pd3==0:
                           sg.Popup('??????q?????????????????????\n???????????????????????????')
                       elif pd3==1:
                           pd4=1
                           sg.Popup('?????????q?????????', font=("times new roman", 15))
                           win5['input'].update(q)
                   elif event5=='??????':
                       sg.Popup('?????????q???????????????????????????????????????\n???????????? \'q=C2H4:1 H2:1 O2:1\'\n??????????????????????????????????????????????????????????????????0000(??????)???????????????????????????\n????????????????????????????????????',title='??????',font=("times new roman", 15))
           elif event3 == '??????':
               win3.close()
               win3_active = False
               win2.UnHide()
               break
           elif event3 == '????????????':
               try:
                   pda()
                   gas1.transport_model = 'Multi'
                   cov = gas1.viscosity
                   win3['out2'].update('??????????????????{}'.format(cov))
               except:
                   sg.Popup('???????????????????????????????????????????????????')
           elif event3 == '??????????????????':
               try:
                   pda()
                   # ??????????????????  ??????????????????????????????????????????????????????????????????????????????????????????????????????
                   massdiff = gas1.mix_diff_coeffs_mass
                   win3['out2'].update('????????????????????????{}'.format(massdiff))
               except:
                   sg.Popup('???????????????????????????????????????????????????')
           elif event3 == '????????????':
               try:
                   pda()
                   # ????????????  ?????????????????????????????????????????????1m???????????????????????????????????????1??????K?????C???,???1???????????????1???????????????????????????????????????????????
                   # ????????????/???????????W/m??k???W/m??K,?????????K?????????????????????
                   heatd = gas1.thermal_conductivity
                   win3['out2'].update('??????????????????{}'.format(heatd))
               except:
                   sg.Popup('???????????????????????????????????????????????????')
           elif event3 == '??????':
               win3['T1'].update(300)
               win3['P1'].update(101300)
               for i in range(lens):
                   win3['fywn{}'.format(i)].update(1)
           elif event3 == '??????':
               win3['T1'].update('')
               win3['P1'].update('')
               for i in range(lens):
                   win3['fywn{}'.format(i)].update('')
           elif event3 == '??????...':
               sg.Popup('About this program', 'Version 3.0', 'PySimpleGUI rocks...','By:Solo dance')
           elif event3 == '??????':
               filename = sg.PopupGetFile('file to open', no_window=True)
               os.system(filename)
           elif event3 == '??????':
               sg.Popup('????????????txt????????????????????????????????????')
               savename = sg.PopupGetFile('file to open', no_window=True)
               if savename!='':
                   with open(savename, 'a') as f:
                         f.write(values3['out1'])
                         f.write(values3['out2'])
           elif event3 == '??????':
               try:
                   pda()
                   lay4 = [[sg.Text('?????????????????????????????????')],
                           [sg.ProgressBar(100, orientation='h', size=(20, 20), key=('progressbar'))]
                           ]
                   win4 = sg.Window('?????????...', default_element_size=(40, 10), grab_anywhere=True,
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
                   #  print("?????????CJ?????????",cj_speed,"m/s")
                   # calculate the CJ speed and find the VN point
                   # find the VN point
                   gas = PostShock_fr(cj_speed, P1, T1, q, mech)
                   denref = gas.density  # ??????
                   Tref = gas.T  # VN??????
                   Pref = gas.P  # VN??????
                   uvn = cj_speed * rho1 / denref
                   gamma = gas.cp / gas.cv  # ?????????
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

                   # ?????????
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
                   # M_CJ  Ma=(r*R/M*T)**0.5  R???????????????????????????R=8.314J/mol???K??????
                   #  r???????????????????????????????????????????????????????????????????????????
                   #  M???????????????????????????
                   #  T??????????????????????????????
                   R = ct.gas_constant
                   a = (gamma * R * T1 / (16 + 1 + 12)) ** 0.5
                   Mcj = cj_speed / a
                   #  print ("Mcj" , Mcj)
                   # heat releas  Q/(R*T0) ????????????
                   Q = (Mcj - 1 / Mcj) * (Mcj - 1 / Mcj) * gamma / 2 / (gamma * gamma - 1)
                   for i in range(70, 100):
                       event4, values4 = win4.read(timeout=10)
                       progress_bar.UpdateBar(i + 1)
                   #  print("Q/RT0", Q)
                   # cov = gas1.viscosity
                   #  print("???????????????",cov)
                   # d=gas.density
                   # cov1=gas1.viscosity
                   # d1=gas1.density
                   #  print("???????????????",heatd)
                   win3['out1'].update(
                       '??????????????????\n  CJ??????:{}\n  Tvn:{}\n Dvn:{}\n  gamma:{}\n  tignition:{}\n  texplosion:{}\n  activation temperature:{}\n  Mcj:{}\n  Q/RTO:{}\n  '.format(
                           cj_speed, Tref, denref, gamma, tignition, texplosion, activationtemperature, Mcj, Q
                       ))
                   win4.close()
               except:
                   sg.Popup('???????????????????????????????????????????????????')


