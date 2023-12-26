from math import log10

_si_prefix = {
     'milimolar':3,
     'micromolar':6,
     'nanomolar':9,
     'picomolar':12
}

def ic50_to_pic50(ic50,mode='nanomolar'):
     try:
        x = _si_prefix[mode.lower()]
        pic50 = x - log10(ic50)
        return pic50
     except KeyError:
         raise ValueError(f'mode must be in {_si_prefix.keys()}')

def pic50_to_ic50(pic50,mode='nanomolar'):
     try:
        x = _si_prefix[mode.lower()]
        ic50 = 10 ** (x - pic50)
        return ic50
     except KeyError:
         raise ValueError(f'mode must be in {_si_prefix.keys()}')

    

    



