import re

def index(i,j,k):

    '''
    Fargo3d index calculator.
    Input: Strings i,j,k, with the value of the desired index on each direction.
    Output: The monodimentional fargo3d index.
    '''

    value = ''
    print i,j,k
    #Trivial option
    if i == 'i' and j == 'j' and k == 'k':
        value += 'l'
        return value
    if i == 'i' or j == 'j' or k == 'k':
        value += 'l'
    
    x = re.match("\w([+-])(\d+)?",i)
    y = re.match("\w([+-])(\d+)?",j)
    z = re.match("\w([+-])(\d+)?",k)

    if x != None:
        if int(x.group(2)) >= 2:
            print '\nError! The allowed displacement in i direction is up to +/- 1\n'
            return
        if x.group(1) == '+':
            value += 'lxp'
        if x.group(1) == '-':
            value += 'lxm'
    
    if y != None:
        if(y.group(2) == '1'):
            value += y.group(1) + 'Nx'
        else:
            value += y.group(1) + y.group(2) + '*Nx'

    if z != None:
        if(z.group(2) == '1'):
            value += z.group(1) + 'Stride'
        else:
            value += z.group(1) + z.group(2) + '*Stride'
    return value
