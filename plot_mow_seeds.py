#!/usr/bin/env python3
import os,sys,re

def fact(n):
    res = 1
    for i in range(2,n+1):
        res *= i
    return res

def binomial(n, k):
    return fact(n) // (fact(k) * fact(n - k))


from fractions import Fraction


# Model "file_4iupac" below 
# - First state 0 (final = 1) is not used ... transitions for the 16 letters are set to 0 (this is not important)
#
# - Second state 1 (final = 0) is self-looping with one transition per letter
#   transition matrix is :
#     A    C    G    T
#  A  3x   xp   xp   xp
#  C  xp   3x   xp   xp
#  G  xp   xp   3x   xp
#  T  xp   xp   xp   3x
#
# but since its integer based computation and the sum is 12 (and not 1), this has to be divided by 12**alignment_length after
# the full computation
#

file_4iupac = """
2   0 1    0 0
           1 0
           2 0
           3 0
           4 0
           5 0
           6 0
           7 0
           8 0
           9 0
           10 0
           11 0
           12 0
           13 0
           14 0
           15 0
    1 0    0 1     1 3 * x
           1 1     1 xp
           2 1     1 xp
           3 1     1 xp
           4 1     1 xp
           5 1     1 3 * x
           6 1     1 xp
           7 1     1 xp
           8 1     1 xp
           9 1     1 xp
           10 1     1 3 * x
           11 1     1 xp
           12 1     1 xp
           13 1     1 xp
           14 1     1 xp
           15 1     1 3 * x
"""

def write_model(filename):
    if not os.path.exists(filename):
        with open(filename,"wt") as f:
            f.write(file_4iupac)


# (1) Extracting seed coeeficients from iedera output
def seeds_coefficients(seed_pattern = "NNNnNnnNnNnnNNnNNN", alignment_length = 64):
    """
    Extracting the iedera polynomial coefficients or counts from iedera (they are only given with the "-pF" option applied)

    This function runs "iedera (>= v1.7)" on the given "seed pattern (iupac notation)" with an "alignment_length", and reads its output.
    It returns a tuple of 2 elements, each element being a list of len = alignment_length + 1 :
    - the   first list gives the    number of alignments with "i" matches (index i) detected by the iupac seed
    - the  second list gives the frequency of "alignments detected by the seed"  /  "alignments in total"

    @param str : the seed pattern being used, as a iedera string on the 27 iupac ("ACGTRrYySsWwKkMmBbDdHhVvn@N") extended symbols for seeds, used by "iedera --iupac")
    @param int : the alignment length
    @return tuple(list<Fraction>,list<Fraction>) : the tuple with two lists
    """
    write_model("file_4iupac.txt")
    err = os.system('./src/iedera -iupac -l '+str(alignment_length)+' -s 1,64 -m "'+seed_pattern+'" -pF file_4iupac.txt -o iedera_output.txt')
    if err != 0:
        sys.exit("./src/iedera not running (please try to run \"./configure; make\" before)")

    file = open('iedera_output.txt','r')
    for line in file:
        matches_pos = re.findall('([0-9]*) \* x \^ ([0-9]*)', line, re.DOTALL)
        # transform tuples of strings into tuples of int
        if matches_pos != []:
            matches_pos_indiced      = [0 for _ in range(alignment_length+1)]
            matches_pos_indiced_frac = [0 for _ in range(alignment_length+1)]
            # transform each tuple by its pos --> value, AND divide by 12**alignment_length, since all the probs per state are summing to 12 ;-)
            for c1,c2 in matches_pos:
                matches_pos_indiced     [int(c2)] = Fraction(int(c1),12**alignment_length)
                matches_pos_indiced_frac[int(c2)] = Fraction(int(c1),12**alignment_length * binomial(alignment_length,int(c2)))

    file.close()
    os.remove('iedera_output.txt')
    return matches_pos_indiced,matches_pos_indiced_frac




# (2) Computation of probabilities on rational numbers
def poly_prob_frac(pos_coef,x):
    """
    Compute the polynomial \sum_i  pos_coef[i] * x ^ i * (1-x) ^ (n-i)
    {where n = len(pos_coef) + 1}
    as a rational number (Fraction)

    @param list<Fraction> : the list of integers coefficients for the polynomial
    @param Fraction       : the value for x as a rational number (Fraction)
    """
    value = Fraction(0)
    n = len(pos_coef)
    for i in range(len(pos_coef)):
       value += pos_coef[i] * (x)**i  *  (1-x)**(n - i - 1)
    return value

# print polynomial
def poly_str(pos_coef):
    """
    Give the the polynomial \sum_i  pos_coef[i] * x ^ i * (1-x) ^ (n-i) as a string

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
def plot_for_seed(plt = None, seed_pattern = "NNNnNnnNnNnnNNnNNN", color_plot = 'gray', alignment_length = 64):
    """
    Plot (and compute before ...) on the current plt the frequency and probability for the seed to hit
    """
    y_pos_coef,y_pos_frac = seeds_coefficients(seed_pattern,alignment_length)

    # (1) Frequency version on rational number values
    x_rational_i1   = [(i1/alignment_length) for i1 in range(alignment_length+1)]
    plt.plot(x_rational_i1, y_pos_frac, marker='.',linestyle = 'None', color=color_plot, label=''+seed_pattern+' frequency')

    # (2) Probabilities on rational numbers values
    x_rational_i2    = [Fraction(i2,delta_range) for i2 in range(delta_range+1)]
    y_rational_poly  = [poly_prob_frac(y_pos_coef,x_rational) for x_rational in x_rational_i2]
    plt.plot(x_rational_i2, y_rational_poly, color=color_plot, label=''+seed_pattern+' probability')

    print(seed_pattern+" : "+poly_str(y_pos_coef))




# MAIN PROGRAMM
import matplotlib.pyplot as plt

alignment_length = 64

# Check arguments
if len(sys.argv) <= 1:
    # Compute the plots for 3 test seeds
    plot_for_seed(plt, seed_pattern = "NNNNNnnNNNnNNnNnnNnNNNnNNNN", color_plot = 'lightgray',  alignment_length = alignment_length) # "NNNNNnnNNNnNNnNnnNnNNNnNNNN" is of weight 19
    plot_for_seed(plt, seed_pattern = "RYNNNNNNNNNNNN",              color_plot = 'darkblue',   alignment_length = alignment_length) # weight 13 but with a  ~ "1/4 x 1/4"  "index x query" sparsity  (so at least ~ "weight 15" considered)
    plot_for_seed(plt, seed_pattern = "RYNNNNnnNNNnnNNNNN",          color_plot = 'darkorange', alignment_length = alignment_length) # (same)
    print("\033[93m This is an example : put your seeds as command line arguments !!! \033[0m")
    print("\033[93m $"+sys.argv[0]+" \"RYNNNNNNNnnnNNNNN\"  \"RYNNNNnnNNNnnNNNNN\" \"RRYNYNNNNNNNNNN\" \033[0m")
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
