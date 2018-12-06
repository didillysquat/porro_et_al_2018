''' This file will be about processing the data that Porro has sent me and the output from the SP submission.
We want to end up with both a redundant fasta that we can feed into the network making and a set of .fasta and .names
files that we can work with. It may also be worth taking into account that Porro may not have taken account of whether
the SP sequences are subsets of her sequences. We'll have to check that but lets do that once we've got everything
else sorted first'''

def pre_processing():
    # this is just to create a fasta for the t50 sequnces
    # let's start by producing a fasta from this .nex file
    nex_input_path = '/Users/humebc/Google_Drive/projects/barbara_forcioli/seqITS_T50_181129.nex'

    # I've already outputed csv abundance dicts so we just need the unique sequencs which seem to be in this .nex
    with open(nex_input_path, 'r') as f:
        nex_file = [line.rstrip() for line in f]

    t50_fasta = []
    seq_start = False
    for i, line in enumerate(nex_file):
        if 'Matrix' in line:
            first_i = i + 1
            seq_start = True

        if ';' in line and seq_start:
            # then i-1 is the last seq line
            last_i = i
            break

    for i in range(first_i, last_i):
        t50_fasta.extend([nex_file[i].split('\t')[0], nex_file[i].split('\t')[1]])

    # here we have the t50_fasta populated
    # now lets write it out including the '>'
    with open('/Users/humebc/Google_Drive/projects/barbara_forcioli/t50_seqs.fasta', 'w') as f:
        for i in range(0, len(t50_fasta), 2):
            f.write('>{}\n'.format(t50_fasta[i]))
            f.write('{}\n'.format(t50_fasta[i + 1].replace('-','')))




def create_names_files():
    # the purpose of this will be to have a fasta that has multiple versions of the same sequence to represent the
    # size of the node for that sequence

    # create a dict from the t50 fasta
    t50_fasta_path = '/Users/humebc/Google_Drive/projects/barbara_forcioli/t50_seqs_corrected.fas'
    with open(t50_fasta_path, 'r') as f:
        t50_file = [line.rstrip() for line in f]
    t50_fasta_dict = {}
    for i in range(0, len(t50_file), 2):
        t50_fasta_dict[t50_file[i][1:]] = t50_file[i+1]

    # create a dict from the sp fasta
    sp_fasta_path = '/Users/humebc/Google_Drive/projects/barbara_forcioli/2018-10-31_16-03-45.679838.DIVs.fasta'
    with open(sp_fasta_path, 'r') as f:
        sp_file = [line.rstrip() for line in f]
    sp_fasta_dict = {}
    for i in range(0, len(sp_file), 2):
        sp_fasta_dict[sp_file[i][1:]] = sp_file[i+1]

    # create the t50 abund dict
    t50_abund_path = '/Users/humebc/Google_Drive/projects/barbara_forcioli/t50_abund.csv'
    t50_abund_dict = {}
    with open(t50_abund_path, 'r') as f:
        t50_abund_file = [line.rstrip() for line in f]
    for line in t50_abund_file:
        t50_abund_dict[str(line.split(',')[0])] = int(line.split(',')[1])

    # create the sp abund dict
    sp_abund_path = '/Users/humebc/Google_Drive/projects/barbara_forcioli/sp_abund.csv'
    with open(sp_abund_path, 'r') as f:
        sp_abund_dict = {str(line.split(',')[0]): int(line.split(',')[1]) for line in f}

    # now create the names file that we will write out and that will be used as input to the network making
    t50_names_file = []
    for seq_name, abund in t50_abund_dict.items():
        count = 0
        # generate a list of sequence that will be joined by commars
        abund_list_items = [seq_name]
        for i in range(abund - 1):
            abund_list_items.append('{}_{}'.format(seq_name, count))
            count += 1
        t50_names_file.append('{}\t{}'.format(seq_name, ','.join(abund_list_items)))

    sp_names_file = []
    for seq_name, abund in sp_abund_dict.items():
        count = 0
        # generate a list of sequence that will be joined by commars
        abund_list_items = [seq_name]
        for i in range(abund - 1):
            abund_list_items.append('{}_{}'.format(seq_name, count))
            count += 1
        sp_names_file.append('{}\t{}'.format(seq_name, ','.join(abund_list_items)))

    # here we have the redundant fasta for each of the t50 and sp created.
    # now write them out
    with open('/Users/humebc/Google_Drive/projects/barbara_forcioli/t50.names', 'w') as f:
        for line in t50_names_file:
            f.write('{}\n'.format(line))

    with open('/Users/humebc/Google_Drive/projects/barbara_forcioli/sp.names', 'w') as f:
        for line in sp_names_file:
            f.write('{}\n'.format(line))

# now we have the redundant fasta, we need to look at how the sequences are related
# something strange is going on as the sums of the t50 seqs add up to 3655 where we end up with like 400 from sp.
def generate_colour_lists():
    # get unique fasta dict for each of the sp and t50
    t50_fasta_path = '/Users/humebc/Google_Drive/projects/barbara_forcioli/t50_seqs_corrected.fas'
    with open(t50_fasta_path, 'r') as f:
        t50_file = [line.rstrip() for line in f]
    t50_fasta_dict = {}
    for i in range(0, len(t50_file), 2):
        t50_fasta_dict[t50_file[i][1:]] = t50_file[i + 1]

    # create a dict from the sp fasta
    sp_fasta_path = '/Users/humebc/Google_Drive/projects/barbara_forcioli/2018-10-31_16-03-45.679838.DIVs.fasta'
    with open(sp_fasta_path, 'r') as f:
        sp_file = [line.rstrip() for line in f]
    sp_fasta_dict = {}
    for i in range(0, len(sp_file), 2):
        sp_fasta_dict[sp_file[i][1:]] = sp_file[i + 1]

    # create the t50 abund dict
    t50_abund_path = '/Users/humebc/Google_Drive/projects/barbara_forcioli/t50_abund.csv'
    t50_abund_dict = {}
    with open(t50_abund_path, 'r') as f:
        t50_abund_file = [line.rstrip() for line in f]
    for line in t50_abund_file:
        t50_abund_dict[str(line.split(',')[0])] = int(line.split(',')[1])

    # create the sp abund dict
    sp_abund_path = '/Users/humebc/Google_Drive/projects/barbara_forcioli/sp_abund.csv'
    with open(sp_abund_path, 'r') as f:
        sp_abund_dict = {str(line.split(',')[0]): int(line.split(',')[1]) for line in f}

    # first lets look at a quick permutation
    t50_match_seqs = []
    sp_match_seqs = set()
    duplicate_match = 0
    for sp_seq_name, sp_seq in sp_fasta_dict.items():
        # for each of the sp seqs, lets see how many of the t50 seqs they match
        for t50_seq_name, t50_seq in t50_fasta_dict.items():
            if sp_seq in t50_seq:
                sp_match_seqs.add(sp_seq_name)
                if t50_seq_name in t50_match_seqs:
                    duplicate_match += 1
                    print('duplicate match {}'.format(t50_seq_name))
                else:
                    t50_match_seqs.append(t50_seq_name)
                print('SP seq {} matches t50 seq {}'.format(sp_seq_name, t50_seq_name))
    print('matches found for {} t50 seqs'.format(len(t50_match_seqs)))
    print('matches found for {} sp seqs'.format(len(sp_match_seqs)))

    # lets write out a colour dictionary so that we know which colours we want the nodes to be
    # lets have it in the format: seq_name,colour_hex\n


    with open('/Users/humebc/Google_Drive/projects/barbara_forcioli/t50_colour.csv', 'w') as f:
        for seq_name in t50_abund_dict.keys():
            if seq_name in t50_match_seqs:
                f.write('{},#A9A9A9\n'.format(seq_name))
            else:
                f.write('{},#ffffff\n'.format(seq_name))

    with open('/Users/humebc/Google_Drive/projects/barbara_forcioli/sp_colour.csv', 'w') as f:
        for seq_name in sp_abund_dict.keys():
            if seq_name in sp_match_seqs:
                f.write('{},#A9A9A9\n'.format(seq_name))
            else:
                f.write('{},#ffffff\n'.format(seq_name))


    apples = '/Users/humebc/Google_Drive/projects/barbara_forcioli/t50_grey_list'

# at this point, we have a redundant fasta for each

generate_colour_lists()
