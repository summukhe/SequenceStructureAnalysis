import random
import numpy as np
import pandas as pd


def get_invalid_char( chr_str ):
    list_def =  '$#^~:><'
    for l in list_def:
        if l not in chr_str:
            return l


def run_length_encode( chr_str ):
    assert isinstance(chr_str , list )
    n = len(chr_str)
    chr_str.append( get_invalid_char(chr_str) )
    rle, count = list(), list()
    counter = 1
    for i in range(n):
        if chr_str[i] == chr_str[i+1]:
            counter += 1
        else:
            rle.append(chr_str[i])
            count.append( counter )
            counter = 1
    return pd.DataFrame( { 'char' : rle, 'count': count } )


def random_sequence_generator( n , chr_set ):
    chr_set = list(set(chr_set))
    array = list()
    for i in range(n):
        idx = random.randint( 0, len(chr_set) - 1 )
        array.append(chr_set[idx])
    return array


def scale_values( values, minimum=0., maximum=0.):
    assert isinstance(values, list)
    values = np.array(values)
    mn, mx = np.min(values) , np.max(values)
    values = (maximum - minimum) * (values - mn)/(mx - mn) + minimum
    return values.tolist()