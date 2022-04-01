# -*- coding: utf-8 -*-
# __author__ = "fengyitong 2018-11-20"

import sys


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

def print_flower_letter(letter):
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

if len(sys.argv) != 2:
    exit('%s \"the words you want to say\"'%sys.argv[0])
trans_string = sys.argv[1].split(' ')
while ' ' in trans_string:
    trans_string.remove(' ')

if trans_string:
    for string in trans_string:
        for i in string:
            if not i in abelt:
                exit('i can not tranform %s' % string)
        else:
            print_flower_letter(string)
else:
    exit('please input what you want to say')
