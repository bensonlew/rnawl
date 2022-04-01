# -*- coding: utf-8 -*-
# __author__ = "fengyitong 2018-11-20"


from __future__ import print_function
import argparse


abelt = dict(
    a = '''
          
     /\\
    /  \\
   / /\ \\
  / ____ \\
 /_/    \_\\
           ''',
    b = '''
  ____
 |  _ \\
 | |_) |
 |  _ <
 | |_) |
 |____/
        ''',
    c = '''
    _____
  / ____|
 | |
 | |
 | |____
  \_____|
         ''',
    d =  '''
  _____
 |  __ \\
 | |  | |
 | |  | |
 | |__| |
 |_____/
         ''',
    e = '''
  ______ 
 |  ____|
 | |__   
 |  __|  
 | |____ 
 |______|
        ''',
    f = '''
  ______
 |  ____|
 | |__
 |  __|
 | |
 |_|
         ''',
    g = '''
    _____
  / ____|
 | |  __
 | | |_ |
 | |__| |
  \_____|
         ''',
    h = '''
  _    _
 | |  | |
 | |__| |
 |  __  |
 | |  | |
 |_|  |_|
         ''',
    i = '''
  _____
 |_   _|
   | |
   | |
  _| |_
 |_____|
         ''',
    j  = '''
       _
      | |
      | |
  _   | |
 | |__| |
  \____/
         ''',
    k = '''
  _  __
 | |/ /
 | ' /
 |  <
 | . \\
 |_|\_\\
       ''',
    l = '''
  _
 | |
 | |
 | |
 | |
 |_|
       ''',
    m = '''
  __  __
 |  \/  |
 | \  / |
 | |\/| |
 | |  | |
 |_|  |_|
         ''',
    n = '''
  _   _
 | \ | |
 |  \| |
 | . ` |
 | |\  |
 |_| \_|
        ''',
    o = '''
   ____
  / __ \\
 | |  | |
 | |  | |
 | |__| |
  \____/
         ''',
    p = '''
  _____
 |  __ \\
 | |__) |
 |  ___/
 | |
 |_|
         ''',
    q = '''
   ____
  / __ \\
 | |  | |
 | |  | |
 | |__| |
  \___\_\\
         ''',
    r = '''
  _____
 |  __ \\
 | |__) |
 |  _  /
 | | \ \\
 |_|  \_\\
         ''',
    s = '''
   _____
  / ____|
 | (___
  \___ \\
  ____) |
 |_____/
         ''',
    t = '''
  _______
 |__   __|
    | |
    | |
    | |
    |_|
          ''',
    u = '''
  _    _
 | |  | |
 | |  | |
 | |  | |
 | |__| |
  \____/
         ''',
    v = '''
 __      __
 \ \    / /
  \ \  / / 
   \ \/ /  
    \  /   
     \/    
           ''',
    w = '''
 __          __
 \ \        / /
  \ \  /\  / /
   \ \/  \/ /
    \  /\  /
     \/  \/
               ''',
    x = '''
 __   __
 \ \ / /
  \ V /
   > <
  / . \\
 /_/ \_\\
        ''',
    y = '''
 __     __
 \ \   / /
  \ \_/ /
   \   /
    | |
    |_|
          ''',
    z = '''
  ______
 |___  /
    / /
   / /
  / /__
 /_____|
        '''
)

bgcolor = {
'black': '40',
'red': '41',
'green': '42',
'yellow': '43',
'blue': '44',
'purple': '45',
'cyan': '46',
'white': '47'
}

fgcolor = {
'black': '30',
'red': '31',
'green': '32',
'yellow': '33',
'blue': '34',
'purple': '35',
'cyan': '36',
'white': '37'
}

model = {
'default': '0',
'highlight': '1',
'underline': '4',
'blink': '5',
'anti': '7',
'hide': '8'
}

def print_flower_letter(letter, bg, fg, mo):
    print('\033[%s;%s;%sm'%(model[mo],fgcolor[fg],bgcolor[bg]),end='\n')
    print_info = [abelt[x] for x in letter.lower()]
    length_line = [len(x.strip().split('\n')) for x in print_info]
    for n in range(max(length_line) + 1):
        info = ' '
        for i in print_info:
            i = i.strip('\n').split('\n')
            try:
                info += i[n].replace('\n', ' ') + ' '*(len(i[-1].replace('\n', ' '))-len(i[n].replace('\n', ' '))) + ' '
            except:
                info += ' '*len(i[-1])
        print(info)
    print('\033[0m',end='')



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="The script flowerlly print your input")
    parser.add_argument("-words", type=str, required=True, help="the words you want to print")
    parser.add_argument("-bgcolor", type=str, default='black', help='the background color, you can choose black,red,green,yellow,blue,purple,cyan,white')
    parser.add_argument("-fgcolor", type=str, default='white', help='the foreground color, you can choose black,red,green,yellow,blue,purple,cyan,white')
    parser.add_argument("-model", type=str, default='default', help="the model of print, you can choose default, highlight, underline,blink,anti,hide")

    args = parser.parse_args()

    trans_string = args.words.split(' ')
    while u'' in trans_string:
        trans_string.remove('')
    
    bg = args.bgcolor
    if bg not in bgcolor:
        bg = 'black'
    
    fg = args.fgcolor
    if fg not in fgcolor:
        bg = 'white'
    
    mo = args.model
    if mo not in model:
        mo = 'default'
    
    
    if trans_string:
        for string in trans_string:
            for i in string:
                if not i in abelt:
                    exit('i can not tranform %s' % string)
            else:
                print_flower_letter(string, bg, fg, mo)
    else:
        exit('please input what you want to say')

