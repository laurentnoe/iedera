#!/usr/bin/env python3
import os,sys,re


# (1) Extracting seed coeeficients from iedera output
def seeds_coefficients(seed_pattern = "###-#--#-#--##-###", alignment_length = 64):
    """
    Extracting the iedera polynomial coefficients or counts from iedera (they are only given with the "-p" option applied)

    This function runs "iedera (>= v1.6)" on the given "seed pattern" with an "alignment_length", and reads its output.
    It returns a tuple of 4 elements, each element being a list of len = alignment_length + 1 :
    - the   first list gives the    number of alignments with "i" matches (index i) detected by the seed
    - the  second list gives the    number of alignments with "i" matches (index i)   missed by the seed
    - the   third list gives the    number of alignments with "i" matches (index i) in total (as the sum of two previous items)
    - the  fourth list gives the frequency of "alignments detected by the seed"  /  "alignments in total"

    @param str : the seed pattern being used (as a iedera string on "#-")
    @param int : the alignment length
    @return tuple(list<int>,list<int>,list<int>,list<float>)  : the tuple with four lists
    """
    list_tuples = ()
    err = os.system('./src/iedera -spaced -l '+str(alignment_length)+' -s 1,64 -m "'+seed_pattern+'" -p -o iedera_output.txt')
    if err != 0:
        sys.exit("./src/iedera not running (please try to run \"./configure; make\" before)")
        
    file = open('iedera_output.txt','r')
    for line in file:
        #print(line,end='')
        matches_pos = re.findall('([0-9]*),1=([0-9]*);', line, re.DOTALL)
        matches_neg = re.findall('([0-9]*),0=([0-9]*);', line, re.DOTALL)
        # transform tuples of strings into tuples of int
        if matches_pos != [] and  matches_neg != []:
            matches_pos_indiced = [0 for _ in range(alignment_length+1)]
            matches_neg_indiced = [0 for _ in range(alignment_length+1)]
            # transform each tuple by its pos --> value
            for c1,c2 in matches_pos:
                matches_pos_indiced[int(c1)] = int(c2)
            for c1,c2 in matches_neg:
                matches_neg_indiced[int(c1)] = int(c2)
    file.close()
    os.remove('iedera_output.txt')
    # compute, from the number of alignments detected/missed by the seed, the total number of alignments that have "i" mismatches
    # { this must give "binomial(i,alignment_length)" }
    # compute also the relative frequencies of success for such alignments
    matches_tot_indiced = [(matches_pos_indiced[i] + matches_neg_indiced[i]) for i in range(alignment_length+1)]
    matches_frq_indiced = [(matches_pos_indiced[i] / matches_tot_indiced[i]) for i in range(alignment_length+1)]
    
    #print(matches_pos_indiced) # number of alignments detected by the seed
    #print(matches_neg_indiced) # number of alignments missed   by the seed
    #print(matches_tot_indiced) # number of alignments in total (binomial coefficient)
    #print(matches_frq_indiced) # frequency of alignments detected / alignments in total
    
    return matches_pos_indiced,matches_neg_indiced,matches_tot_indiced,matches_frq_indiced




# (2) Computation of probabilities performed by including smallest values first
# basic float version
def poly_prob_float(pos_coef,x):
    """
    Compute the polynomial \\sum_i  pos_coef[i] * x ^ i * (1-x) ^ (n-i)
    {where n = len(pos_coef) + 1}
    as a floating point value


    @param list<int> : the list of integers coefficients for the polynomial
    @param float : the value for x as a floating point number
    """
    value = 0
    n = len(pos_coef)
    for i in range(len(pos_coef)):
       value += pos_coef[i] * (x)** i  *  (1-x)**(n - i - 1) 
    return value


# rational version
from fractions import Fraction
def poly_prob_frac(pos_coef,x):
    """
    Compute the polynomial \\sum_i  pos_coef[i] * x ^ i * (1-x) ^ (n-i)
    {where n = len(pos_coef) + 1}
    as a rational number (Fraction)

    @param list<int> : the list of integers coefficients for the polynomial
    @param float : the value for x as a rational number (Fraction)
    """
    value = Fraction(0)
    n = len(pos_coef)
    for i in range(len(pos_coef)):
       value += Fraction(pos_coef[i]) * (x)**i  *  (1-x)**(n - i - 1)
    return value

# str output polynomial
def poly_str(pos_coef):
    """
    Give the the polynomial \\sum_i  pos_coef[i] * x ^ i * (1-x) ^ (n-i) as a string

    @param list<int> : the list of integers coefficients for the polynomial
    @return str      : the representation of this polynomial
    """
    value = ""
    n = len(pos_coef)
    for i in range(len(pos_coef)):
       value += str(pos_coef[i]) + " * x^" + str(i) + " * (1-x)^" + str(n - i - 1)
       if i < len(pos_coef)-1:
           value += " + "
    return value


# (3) plotting scripts for one seed
delta_range = 1000
def plot_for_seed(plt = None, seed_pattern = "###-#--#-#--##-###", color_plot = 'gray', alignment_length = 64):
    """
    Plot (and compute before ...) on the current plt the frequency and probability for the seed to hit    
    """
    # (1) frequency version
    x_i1 = [(i1/alignment_length) for i1 in range(alignment_length+1)]
    y_pos,y_neg,y_tot,y_freq = seeds_coefficients(seed_pattern,alignment_length)
    plt.plot(x_i1, y_freq, marker='.',linestyle = 'None', color=color_plot, label=''+seed_pattern+' frequency')

    # (2.1) Probabilities on floating-point values (precise ?)
    x_float_i2 = [i/delta_range for i in range(delta_range+1)]
    y_float_i2 = [poly_prob_float(y_pos,x_float) for x_float in x_float_i2]
    plt.plot(x_float_i2,y_float_i2, color=color_plot, label=''+seed_pattern+' probability')

    # (2.2) Probabilities on rational numbers values (precise !)
    ## x_rational_i2 = [Fraction(i,delta_range) for i in range(delta_range+1)]
    ## y_rational_i2 = [poly_prob_frac(y_pos,x_rational) for x_rational in x_rational_i2]
    ## plt.plot(x_rational_i2, y_rational_i2, color=color_plot, label=''+seed_pattern+' probability (rational)')

    print(seed_pattern+" : "+poly_str(y_pos))




# MAIN PROGRAMM
import matplotlib.pyplot as plt

alignment_length = 64

# Check arguments
if len(sys.argv) <= 1:
    # Compute the plots for two test seeds
    plot_for_seed(plt, seed_pattern = "###########",        color_plot = 'darkblue',   alignment_length = alignment_length)
    plot_for_seed(plt, seed_pattern = "###-#--#-#--##-###", color_plot = 'darkorange', alignment_length = alignment_length)
    print("\033[93m This is an example : put your seeds as command line arguments !!! \033[0m")
    print("\033[93m $"+sys.argv[0]+" \"############\"  \"###-#--###-#--###-#\" \"####-#-##--####-#-##,#-##--####-#-##--####\" \033[0m")
else:
    # Compute the plots for several given seeds
    # (https://matplotlib.org/stable/gallery/color/named_colors.html)
    colors_plot = ['darkblue','darkorange','gold','orangered','grey','seagreen']
    for i in range(1,len(sys.argv)):
        plot_for_seed(plt, seed_pattern = sys.argv[i], color_plot = colors_plot[(i-1) % len(colors_plot)], alignment_length = alignment_length)

# Plotting all that stuff
plt.xlim([0,1])
plt.xlabel("alignment match frequency/probability")
plt.ylim([-0.01,+1.01])
plt.ylabel("1st hit frequency/probability")

plt.title("1st hit frequency/probability (alignment length = "+str(alignment_length)+")")
plt.legend(loc="upper left")

plt.show()
