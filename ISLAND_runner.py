from Bio import SearchIO
import pandas as pd
import time
import os
import itertools
import numpy
import subprocess
from Bio import SeqIO
from Bio import SeqRecord
from uuid import uuid1
from multiprocessing import Pool

def chunker(xs, n):

    '''Split the list, xs, into n evenly sized chunks'''
    L = len(xs)
    assert 0 < n <= L
    s, r = divmod(L, n)
    t = s + 1
    return ([xs[p:p+t] for p in range(0, r*t, t)] +
            [xs[p:p+s] for p in range(r*t, L, s)])

def location_in_pair(location, pair, threshold):

    lower = min(pair)
    upper = max(pair)

    if upper > location > lower:

        return True

    elif location < lower:

        if abs(upper - location) <= threshold:

            return True

        else:

            return False

    elif location > upper:

        if abs(lower - location) <= threshold:

            return True

        else:

            return False

def logic_brute_force(dictionary, mand, accessory, total_needed, window):

    if all([name in dictionary.keys() for name in mand]) and len(mand) > 0:

        mand_gene_list = []
        accessory_gene_list = []

        for gene in mand:
            mand_gene_list.append(dictionary[gene])

        for gene in accessory:
            if gene in dictionary:
                accessory_gene_list.append(dictionary[gene])

        passed_mand_pairs = pairs_shared_between_all_lists(mand_gene_list, window)

        if len(passed_mand_pairs) > 0:

            if len(accessory) > 0:

                for gene_list in accessory_gene_list:

                    for location in gene_list:

                        for pair in passed_mand_pairs:

                            if location_in_pair(location, pair, window):

                                pair.append(location)

                out_list = []

                for pair in passed_mand_pairs:

                    if len(pair) >= total_needed:

                        out_list.append(pair)

                return out_list
            else:

                return passed_mand_pairs

        else:

            return []

    elif len(mand) == 0:

        accessory_gene_list = []

        for gene in accessory:
            if gene in dictionary:
                accessory_gene_list.append(dictionary[gene])

        if accessory_gene_list == []:

            return []

        longest = max(accessory_gene_list, key=len)

        passed_mand_pairs = [[thing] for thing in longest]

        if len(passed_mand_pairs) > 0:

            if len(accessory) > 0:

                for gene_list in accessory_gene_list:

                    for location in gene_list:

                        for pair in passed_mand_pairs:

                            if location_in_pair(location, pair, window):

                                pair.append(location)

                out_list = []

                for pair in passed_mand_pairs:

                    if len(pair) >= total_needed:

                        out_list.append(pair)

                return out_list

    else:

        return []


def find_pairs_two_lists(list1, list2, threshold):

    list1 = sorted(list1)
    list2 = sorted(list2)

    pairs = []

    for i in range(len(list1)):

        for j in range(len(list2)):

            if abs(list1[i] - list2[j]) <= threshold:

                pairs.append([list1[i], list2[j]])

            elif list2[j] > list1[i]:

                break

    return pairs

def list_and_list_of_pairs(input_list, list_of_pairs, threshold):

    input_list = sorted(input_list)

    if len(list_of_pairs) == 0:

        return []

    if isinstance(list_of_pairs[0], int):

        list_of_pairs = [[pair] for pair in list_of_pairs]

    new_pairs = []

    for pair in list_of_pairs:

        for number in input_list:

            if all([abs(number - value) <= threshold for value in pair]) and number not in pair:

                new_pairs.append(pair + [number])

            elif number - min(pair) > threshold:

                break

    return new_pairs

def pairs_shared_between_all_lists(list_of_lists, threshold):

    if len(list_of_lists) == 0:

        return []

    if len(list_of_lists) == 1:

        return [[thing] for thing in list_of_lists[0]]

    out_list = []
    current_list = []

    for item in list_of_lists:

        working_list = []

        for thing in item:

            working_list.append(thing)

        current_list.append(working_list)

    initial_pairs_list = find_pairs_two_lists(current_list[0], current_list[1], threshold)

    running_pairs_list = initial_pairs_list

    for thing in current_list[2:]:
        running_pairs_list = list_and_list_of_pairs(thing, running_pairs_list, threshold)

    out_list += running_pairs_list

    return out_list


class System:

    def __init__(self, name, mand, accessory, number_needed):

        self.name = name
        self.mand = mand
        self.accessory = accessory
        self.number_needed = number_needed


class Gene:

    def __init__(self, name, models, how_to_find):

        self.name = name
        self.models = models
        self.how_to_find = how_to_find

def ISLAND_runner(gembase_file_and_out_dir, phamlite_tmp):
    temp_dir = phamlite_tmp + str(uuid1()) + '/'
    gembase_file = gembase_file_and_out_dir.split(':')[0]
    cores = gembase_file_and_out_dir.split(':')[1]
    out_csv_link = phamlite_tmp + str(uuid1()) + '.csv'
    subprocess.call('mkdir ' + temp_dir, shell=True)
    # reading in the systems
    start = time.time()
    profiles_dir = '/Users/matt/Downloads/converted_ISLAND/'
    profiles = os.listdir(profiles_dir)
    systems_guide = pd.read_csv('/hdd-roo/ericl/ISLAND/ISLAND_updated_master.csv')
    systems_dict = {}
    model_dict = {'Profile': [], 'Search': []}
    for i in range(len(systems_guide['System'])):
        system = systems_guide['System'][i]
        current_gene = systems_guide['Gene'][i]
        find_way = systems_guide['Howto'][i]
        requirement = systems_guide['Requirement'][i]
        needed = systems_guide['Total Needed'][i]
        if system not in systems_dict:
            systems_dict[system] = {'Mandatory': [], 'Accessory': [], 'Total Needed': needed}
        models = []
        for profile in profiles:
            if current_gene in profile:
                model_dict[find_way].append(profiles_dir + profile)
        current_gene_object = Gene(current_gene, models, find_way)
        systems_dict[system][requirement].append(current_gene_object)
    systems = []
    for system in systems_dict:
        systems.append(System(system, systems_dict[system]['Mandatory'], systems_dict[system]['Accessory'],
                              systems_dict[system]['Total Needed']))
    # make the large blast fasta/hmm files for blast and hmmsearch
    blast_maker_command = 'cat ' + profiles_dir + '*.faa > ' + temp_dir + 'blast_proteins.faa'
    subprocess.call(blast_maker_command, shell=True)
    hmm_mega_maker_command = 'cat ' + profiles_dir + '*.hmm > ' + temp_dir + 'hmms_all.hmm'
    subprocess.call(hmm_mega_maker_command, shell=True)
    # cluster the gembase file
    cluster_command = '/usr/local/bin/mmseqs easy-linclust ' + gembase_file + ' ' + temp_dir + 'gembase_cluster ' + temp_dir + ' --min-seq-id 0.99 --threads ' + cores
    subprocess.call(cluster_command, shell=True, stdout=subprocess.DEVNULL)
    # make dict with genome and index and sequence
    current_genome = ''
    current_contig = ''
    index = 1
    position_dict = {}
    contig_genome_dict = {}
    handle = open(gembase_file, 'r')
    for record in SeqIO.parse(handle, 'fasta'):
        try:
            protein = record.id.split('q')[1].split('_')[1]
            if record.id.split('q')[1].split('_')[0] != current_contig:
                current_contig = record.id.split('q')[1].split('_')[0]
                index = 1
            position_dict[protein] = index
            genome = record.id.split('q')[0]
            if current_contig not in contig_genome_dict:
                contig_genome_dict[current_contig] = genome
            index += 1
        except IndexError:
            pass
    handle.close()
    # make the cluster results dict

    cluster_dict = {}

    cluster_data = pd.read_csv(temp_dir + 'gembase_cluster_cluster.tsv', sep='\t', names=['Cluster', 'Member'])

    for i in range(len(cluster_data['Cluster'])):

        cluster = cluster_data['Cluster'][i]
        member = cluster_data['Member'][i]

        if cluster not in cluster_dict:
            cluster_dict[cluster] = []

        cluster_dict[cluster].append(member)


    file_to_search = temp_dir + 'gembase_cluster_rep_seq.fasta'

    # do the hmmsearch

    hmmsearch_command = 'hmmsearch --domtblout ' + temp_dir + 'hmmsearch_out.txt --cpu ' + cores + ' ' + temp_dir + 'hmms_all.hmm ' + file_to_search

    subprocess.call(hmmsearch_command, shell=True, stdout=subprocess.DEVNULL)

    # do the blast search

    blast_db_command = 'makeblastdb -in ' + file_to_search + ' -out ' + temp_dir + 'blast_db -parse_seqids -dbtype prot'
    subprocess.call(blast_db_command, shell=True, stdout=subprocess.DEVNULL)
    blast_search_command = 'blastp -query ' + temp_dir + 'blast_proteins.faa ' + ' -db ' + temp_dir + 'blast_db -out ' + temp_dir + 'blastresults.xml -evalue 0.001 -max_target_seqs 10000000 -outfmt 5 -num_threads ' + cores
    subprocess.call(blast_search_command, shell=True, stdout=subprocess.DEVNULL)


    # parse through hmm first to make genome dict then blast

    contig_dict = {}

    handle = open(temp_dir + 'hmmsearch_out.txt', 'r')

    try:

        for record in SearchIO.parse(handle, 'hmmsearch3-domtab'):

            gene = str(record.id)[0:4]
            model_len = record.seq_len

            for hit in record.hits:

                for hsp in hit.hsps:

                    if hsp.evalue < 0.001 and hsp.query_span / model_len > 0.5:

                        contigs = []

                        contigs.append(hit.id)

                        contigs += cluster_dict[hit.id]

                        for cluster_hit in contigs:

                            contig = cluster_hit.split('q')[1].split('_')[0]
                            protein = cluster_hit.split('q')[1].split('_')[1]

                            if contig not in contig_dict.keys():

                                contig_dict[contig] = {gene: []}
                                contig_dict[contig][gene].append(protein)

                            elif gene not in contig_dict[contig].keys():

                                contig_dict[contig][gene] = [protein]

                            else:

                                contig_dict[contig][gene].append(protein)

        handle.close()
        handle = open(temp_dir + 'blastresults.xml', 'r')

        for record in SearchIO.parse(handle, 'blast-xml'):

            gene = record.id.split(':')[0]
            query_len = record.seq_len

            for hit in record.hits:

                if max([hsp.query_span for hsp in hit.hsps]) / query_len > 0.5:

                    contigs = []

                    contigs.append(hit.id)

                    contigs += cluster_dict[hit.id]

                    for cluster_hit in contigs:

                        contig = cluster_hit.split('q')[1].split('_')[0]
                        protein = cluster_hit.split('q')[1].split('_')[1]

                        if contig not in contig_dict.keys():

                            contig_dict[contig] = {gene: []}
                            contig_dict[contig][gene].append(protein)

                        elif gene not in contig_dict[contig].keys():

                            contig_dict[contig][gene] = [protein]

                        else:

                            contig_dict[contig][gene].append(protein)

        handle.close()
        big_out_dict = {'Contig': []}
        protein_id_contig_out_dict = {'Contig': []}

        for system in systems:
            big_out_dict[system.name] = []
            protein_id_contig_out_dict[system.name] = []

        for contig in contig_dict:

            current_dict = {}
            decoder = {}

            for key in contig_dict[contig]:

                current_dict[key] = []

                for thing in contig_dict[contig][key]:
                    current_dict[key].append(position_dict[thing])
                    decoder[position_dict[thing]] = thing

            for system in systems:

                mand_list = []
                accessory_list = []

                for gene in system.mand:
                    mand_list.append(gene.name)

                for gene in system.accessory:
                    accessory_list.append(gene.name)

                x = logic_brute_force(current_dict, mand_list, accessory_list, system.number_needed, 10)
                big_out_dict[system.name].append(int(len(x) > 0))

                if len(x) > 0:

                    protein_id_contig_out_dict[system.name].append(
                        '&'.join(list(set(':'.join([decoder[thing] for thing in system]) for system in x))))

                else:

                    protein_id_contig_out_dict[system.name].append('')

            big_out_dict['Contig'].append(contig)
            protein_id_contig_out_dict['Contig'].append(contig)

        semi_final_binary_dict = {}
        semi_final_ids_dict = {}

        for i in range(len(big_out_dict['Contig'])):

            contig = big_out_dict['Contig'][i]
            genome = contig_genome_dict[contig]

            if genome not in semi_final_binary_dict:
                semi_final_binary_dict[genome] = {}

            for system in systems:

                if system.name not in semi_final_binary_dict[genome]:
                    semi_final_binary_dict[genome][system.name] = 0

                semi_final_binary_dict[genome][system.name] += big_out_dict[system.name][i]

        for i in range(len(protein_id_contig_out_dict['Contig'])):

            contig = protein_id_contig_out_dict['Contig'][i]
            genome = contig_genome_dict[contig]

            if genome not in semi_final_ids_dict:
                semi_final_ids_dict[genome] = {}

            for system in systems:

                if system.name not in semi_final_ids_dict[genome]:
                    semi_final_ids_dict[genome][system.name] = ''

                semi_final_ids_dict[genome][system.name] += protein_id_contig_out_dict[system.name][i]

        final_out_dict_binary = {'Genome': []}

        for system in systems:
            final_out_dict_binary[system.name] = []

        for genome in semi_final_binary_dict:

            final_out_dict_binary['Genome'].append(genome)

            for system in semi_final_binary_dict[genome]:
                final_out_dict_binary[system].append(semi_final_binary_dict[genome][system])

        final_out_dict_id = {'Genome': []}

        for system in systems:
            final_out_dict_id[system.name] = []

        for genome in semi_final_ids_dict:

            final_out_dict_id['Genome'].append(genome)

            for system in semi_final_ids_dict[genome]:
                final_out_dict_id[system].append(semi_final_ids_dict[genome][system])

        x = pd.DataFrame.from_dict(final_out_dict_id)
        y = pd.DataFrame.from_dict(final_out_dict_binary)
        x.to_csv(out_csv_link.replace('.csv', '_proteins.csv'))
        y.to_csv(out_csv_link)

        print('All done', time.time() - start)

        subprocess.call('rm -r ' +temp_dir, shell=True)

    except IndexError:

        pass

gembase_file_directory = '/hdd-roo/genomes_master/contig_gembases/gpff/'

gembases = [gembase_file_directory + thing + ':8' for thing in os.listdir(gembase_file_directory)]

def run(tmp):
    for thing in gembases:
        p = Pool(12)
        p.map(ISLAND_runner, [thing, tmp])


gembase_dir = '/hdd-roo/genomes_master/contig_gembases/gpff/'
# main_output_dir = '/hdd-roo/ericl/ISLAND/erin_cas_macsy_out_non_clinical/'
#
# subprocess.call('mkdir ' + main_output_dir, shell=True)
#
# for file in os.listdir(gembase_dir):
#
#     if file[-4:] == '.txt':
#
#         macsy_command = 'macsyfinder -w 32 --sequence-db ' + gembase_dir + file + ' --db-type gembase all ' + ' -d /hdd-roo/ericl/CRISPRCasFinder/CasFinder-2.0.3/DEF-SubTyping-2.0.3/ -p /hdd-roo/ericl/CRISPRCasFinder/CasFinder-2.0.3/CASprofiles-2.0.3/ -o ' + main_output_dir + file.replace('.txt', '') + '/'
#         subprocess.call(macsy_command, shell=True)

directory_of_interest = '/hdd-roo/ericl/ISLAND/ISLAND_parts/'

pds_to_cat_protein = []
pds_to_cat_binary = []

for file in os.listdir(directory_of_interest):

    if '_proteins' in file:

        pds_to_cat_protein.append(pd.read_csv(directory_of_interest + file))

big_out_protein = pd.concat(pds_to_cat_protein)


big_out_protein.to_csv('/hdd-roo/ericl/ISLAND/ISLAND_proteins_2_16_redo.csv')
