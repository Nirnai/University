# -*- coding: utf-8 -*-
"""
Created on Wed Mar 02 12:21:41 2016

@author: Nirnai
"""
import numpy as np
from copy import copy

### ToDo: Add taps for more m-sequence orders

# see http://www.newwaveinstruments.com/resources/articles/m_sequence_linear_feedback_shift_register_lfsr.htm#Table%20of%20M-Sequence%20Feedback%20Taps for reference

taps = [

    ### degree 1

    [[]],

    ### degree 2

    [[]],

    ### degree 3

    [

        # 3 stages, 2 taps:  (1 set)

        [3, 2]

    ],

    ### degree 4

    [

        # 4 stages, 2 taps:  (1 set)

        [4, 3],

    ],

    ### degree 5

    [[]],

    ### degree 6

    [[]],

    ### degree 7

    [[]],

    ### degree 8

    [[]],

    ### degree 9

    [[]],

    ### degree 10

    [

        # 10 stages, 2 taps:  (1 set)

        [10, 7],

        # 10 stages, 4 taps:  (10 sets)

        [10, 9, 8, 5],

        [10, 9, 7, 6],

        [10, 9, 7, 3],

        [10, 9, 6, 1],

        [10, 9, 5, 2],

        [10, 9, 4, 2],

        [10, 8, 7, 5],

        [10, 8, 7, 2],

        [10, 8, 5, 4],

        [10, 8, 4, 3],

        # 10 stages, 6 taps:  (14 sets)

        [10, 9, 8, 7, 5, 4],

        [10, 9, 8, 7, 4, 1],

        [10, 9, 8, 7, 3, 2],

        [10, 9, 8, 6, 5, 1],

        [10, 9, 8, 6, 4, 3],

        [10, 9, 8, 6, 4, 2],

        [10, 9, 8, 6, 3, 2],

        [10, 9, 8, 6, 2, 1],

        [10, 9, 8, 5, 4, 3],

        [10, 9, 8, 4, 3, 2],

        [10, 9, 7, 6, 4, 1],

        [10, 9, 7, 5, 4, 2],

        [10, 9, 6, 5, 4, 3],

        [10, 8, 7, 6, 5, 2],

        # 10 stages, 8 taps:  (5 sets)

        [10, 9, 8, 7, 6, 5, 4, 3],

        [10, 9, 8, 7, 6, 5, 4, 1],

        [10, 9, 8, 7, 6, 4, 3, 1],

        [10, 9, 8, 6, 5, 4, 3, 2],

        [10, 9, 7, 6, 5, 4, 3, 2],

    ],

]











def get_taps(degree=10, number=0):

    r"""
    Read taps for M-sequence generation
    Parameters
    ----------
    degree: int Degree of sequence
    number: Choose number of M-sequence
    Returns
    -------
    array_like
        Taps for LFSR
    """
    
    if number < len(taps[degree-1]):

        return np.array(taps[degree-1][number])

    elif number < 2*len(taps[degree-1]):

        number-= len(taps[degree-1])

        ret = degree-np.array(taps[degree-1][number])

        ret[0] += degree

        return ret

    else:

        return None










def m_sequence(m, taps = np.array([]), register = np.array([]), delay=0, bi=True):
    
    r"""
    Generate M-sequence
    Parameters
    ----------
    m: int degree of sequence
    taps: array_like Taps of LFSR (see function get_taps())
    register: array_like Initial state of LFSR (default [1,0,0,...])
    delay: int Delay Sequence (default: 0)
    bi: bool Bipolar output (default: True)
    Returns
    -------
    array_like M-sequence
    """
    
    if not type(taps) == 'numpy.ndarray':

        taps = np.array(taps)

    if not type(register) == 'numpy.ndarray':

        np.array(register)

    if taps.shape == (0,):

        taps = get_taps(m)

    if register.shape == (0,):

        register = np.array([1])

    sequenceLength  =   2**m-1

    seq             =   np.zeros(sequenceLength)

    if register.shape != (m,):
        
        register = copy(register)
        register.resize((m,))
        
    registerStates  =   np.zeros((sequenceLength, m))
    
    for i in xrange(0,sequenceLength):
        
        registerStates[i,:] = register
        # Output = last value of current register
        seq[i]          =   register[m-1]
        # Calculate the feedback value
        feedback        =   np.sum(register[np.array(taps)-1]) % 2
        # Update the register
        register        =   np.concatenate((np.array([feedback]), register[0:-1]))
    if bi:
        seq = -2 * seq + 1
    if delay != 0:
        seq = np.roll(seq, delay)
        
    return seq