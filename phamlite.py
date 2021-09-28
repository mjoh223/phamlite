import dash
import dash_bootstrap_components as dbc
from dash import dcc
from dash import html
from dash.dependencies import Input, Output, State
import pandas as pd
from Bio import SeqIO, SeqUtils
from Bio.Seq import Seq
#from seeker import SeekerFasta
from Bio.SeqRecord import SeqRecord
import os
from itertools import combinations
from Bio.Blast import NCBIXML
import glob
from Bio.Blast.Applications import NcbiblastnCommandline
import shutil
import jsonpickle
import subprocess
import csv
import plotly.graph_objects as go
import random
import string
import base64
import datetime
import io
import json
import numpy as np
from colour import Color

#binb = '/home/jbd_apps/dash_app/app_env/bin'
#external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = dash.Dash(__name__)#, external_stylesheets=external_stylesheets)
server = app.server
def makedir(f):
    directories = ['tmp','faa','fna','blast_out','cluster_data','cluster_out']
    for basename in directories:
        os.mkdir( os.path.join(f, basename) )

class Locus:
    def __init__(self, ann, fna, filename):
        self.orfs = []
        self.tRNAs = []
        self.repeat_regions = []
        self.annotations = ann
        self.fna = fna
        self.filename = filename
        self.phams = []
        self.gc_skew = []
        self.element_type = self.annotations['taxonomy'][0].lower()
        self.accessions = self.annotations['accessions'][0]
    def load_gc_trace(self, z, window):
        trace_list = []
        self.gc_skew = SeqUtils.GC_skew(self.fna, window)
        length = int(len(self.fna))
        bin_size = int(length/window)
        red = Color("red")
        blue = Color("blue")
        color_scale = list(red.range_to(blue, 1000))
        for i, b in enumerate( range(0, length, bin_size) ):
            if self.gc_skew[i] > .5:
                color = 'blue'
            else:
                color= 'red'
            x, y = draw_shape(b, bin_size, 'nostrand' ,z, 0.2, 'na', 'na')
            trace = go.Scatter(x=x, y=y, marker=dict(size=1), opacity=1,fill='toself', fillcolor=color, line_color=color )
            trace_list.append(trace)
        return trace_list

    def toJSON(self):
        return json.dumps(self, default=lambda o: o.__dict__,
            sort_keys=True, indent=4)
    def add_feature(self,feature):
        if feature.type == 'CDS':
            self.orfs.append(feature)
        if feature.type == 'tRNA':
            print(feature)
            self.tRNAs.append(feature)
        if feature.type == 'repeat_region':
            self.repeat_regions.append(feature)
    def getPhams(self):
        pham_set = set()
        for orf in self.orfs:
            start, stop, strand = orf.location.start, orf.location.end, orf.location.strand
            id = '|'.join([self.accessions, orf.id[0]]) #orf.id[0]
            loc = np.where(pham_df.member == id)[0]
            pham = pham_df.iloc[loc,0]
            pham_set.add(pham.to_string(index=False))
        self.phams = pham_set
        return pham_set

    def load_orf_trace(self, pham_color_dict, pham_df, z=0, h=0.2):
        self.firstorf = np.min([x.location.start for x in self.orfs])
        self.lastorf = np.max([x.location.start for x in self.orfs])
        fill='toself'
        trace_list = []
        for orf in self.orfs:
            start, stop, strand = orf.location.start, orf.location.end, orf.location.strand
            id = '|'.join([self.accessions, orf.id[0]])
            loc = np.where(pham_df.member == id)[0]
            color = pham_color_dict[list(pham_df.iloc[loc,0])[0]]
            linecolor = 'black'
            opacity=1
            # if orf.getProduct().lower() != 'hypothetical protein':
            #     color = "gray"
            #     opacity = 0.3
            x, y = draw_shape(start, stop, strand, z, h, self.firstorf, self.lastorf)
            trace = go.Scatter(x=x, y=y, marker=dict(size=1), opacity=opacity,fill=fill, fillcolor=color, line_color=linecolor, text='{}|{}'.format(orf.getProduct(), self.annotations['organism']),hoverinfo='text' )
            trace_list.append(trace)
        return trace_list
    def load_trna_trace(self,z=0,h=0.2):
        self.firstorf = np.min([x.location.start for x in self.orfs])
        self.lastorf = np.max([x.location.start for x in self.orfs])
        trace_list = []
        if len(self.tRNAs) > 0:
            for trna in self.tRNAs:
                print(trna)
                print(dir(trna.location))
                if trna.location is None:
                    continue
                start, stop, strand = trna.location.start, trna.location.end, trna.location.strand
                linecolor = 'gold'
                opacity = 1
                x, y = draw_trna(start,stop,strand,z,h,self.firstorf,self.lastorf)
                trace = go.Scatter(x=x, y=y, marker=dict(size=1), opacity=opacity, fill='toself', fillcolor='gold', line_color='gold', text='{}|{}'.format(''.join(trna.product), self.annotations['organism']),hoverinfo='text')
                trace_list.append(trace)
        return trace_list
    def load_repeat_trace(self,z=0,h=0.2):
        self.firstorf = np.min([x.location.start for x in self.orfs])
        self.lastorf = np.max([x.location.start for x in self.orfs])
        trace_list = []
        if len(self.repeat_regions) > 0:
            for rp in self.repeat_regions:
                start, stop, strand = rp.location.start, rp.location.end, rp.location.strand
                linecolor = 'pink'
                opacity = 1
                x, y = draw_repeat(start,stop,strand,z,h,self.firstorf,self.lastorf)
                trace = go.Scatter(x=x, y=y, marker=dict(size=1), opacity=opacity, fill='toself', fillcolor='pink', line_color='pink', text='{}|{}'.format(''.join(rp.rpt_family), self.annotations['organism']),hoverinfo='text')
                trace_list.append(trace)
        return trace_list

    def load_syn_trace(self, blast_di, order, phages, current_h=0):
        shade_trace_list, boundary_left_list, boundary_right_list = [], [], []
        for comparison, matches in blast_di.items():
            if len(matches) > 0 and self.accessions in comparison:
                for match in matches.itertuples():
                    try:
                        target_h = [x[1] for x in order if x[0] in match.subject][0]
                        source_h = [x[1] for x in order if x[0] in match.query][0]
                    except:
                        break
                    if abs(target_h - source_h) > 1:
                        break
                    source, target = match.query, match.subject #am i the source or the target?
                    if self.accessions in source:
                        self.whoami = 'source'
                    else:
                        self.whoami = 'target'
                    source_start, source_end, target_start, target_end, percent_id = match.q_start,match.q_end,match.s_start,match.s_end,match.identity
                    shade = 'purple'
                    if percent_id > 90:
                        shade = 'green'
                    end = len(self.fna)
                    #if source_h in [0]:
                    #    accession = order[source_h][0]
                    #    end = [len(x.fna) for x in phages if x.accessions == accession][0]
                    #    source_start = (end-source_start)
                    #    source_end = (end-source_end)
                    #if target_h in [0]:
                    #    accession = order[target_h][0]
                    #    end = [len(x.fna) for x in phages if x.accessions == accession][0]
                    #    target_start = (end-target_start)
                    #    target_end = (end-target_end)
                    x=(source_start, target_start, target_end, source_end, source_start)
                    y=(source_h, target_h, target_h, source_h, source_h)

                    shade_trace = go.Scatter(x=x,y=y,marker=dict(size=1),fill='toself',fillcolor=shade,line_color=shade,opacity=.1,text='{}%'.format(percent_id),hoverinfo='text')
                    shade_trace_list.append(shade_trace)

                    x=(source_start, target_start, target_end, source_end, source_start)
                    y=(source_h, target_h, target_h, source_h, source_h)
                    shade_trace = go.Scatter(x=x,y=y,marker=dict(size=1),fill='toself',fillcolor=shade,line_color=shade,opacity=.1,text='{}%'.format(percent_id),hoverinfo='text')
                    #boundary_left = go.Scatter(x=(source_start, target_start),y=(source_h, target_h),line_color = shade, mode = 'lines',opacity=.5)
                    #boundary_right = go.Scatter(x=(target_end, source_end),y=(target_h, source_h),line_color = shade, mode = 'lines',opacity=.5)
                    shade_trace_list.append(shade_trace)
                    # boundary_left_list.append(boundary_left)
                    # boundary_right_list.append(boundary_right)


        return shade_trace_list#, boundary_left_list, boundary_right_list

class TRNA(object):
    def __init__(self, feature):
        self.type = feature.type
        self.location = feature.location
        self.qualifiers = feature.qualifiers
        try:
            self.product = feature.qualifiers['product']
        except:
            self.product = 'unknown tRNA'
        try:
            self.note = feature.qualifiers['note']
        except:
            self.note = 'unknown ncRNA or tRNA'

class RepeatRegion(object):
    def __init__(self, feature):
        self.type = feature.type
        self.location = feature.location
        self.qualifiers = feature.qualifiers
        try:
            self.rpt_family = feature.qualifiers['rpt_family']
        except:
            self.rpt_family = 'unknown repeat'
        try:
            self.rpt_unit_seq = feature.qualifiers['rpt_unit_seq']
        except:
            self.rpt_unit_seq = 'unknown repeat sequence'

class Orf(object):
    def __init__(self, feature,i):
        self.type = feature.type
        self.location = feature.location
        self.qualifiers = feature.qualifiers
        try:
            self.id = self.qualifiers['protein_id']
        except:
            self.id = 'protein_{}'.format(i)
    def getSeq(self):
        try:
            protein = self.qualifiers['translation']
        except:
            print(self.id)
            protein = ['']
        return protein
    def getProduct(self):
        return self.qualifiers.get('product', ['unknown'])[0]
    def getPham(self, pham_df):
        phams = set(pham_df.representative)
        rgb_values = ['rgb{}'.format(tuple(np.random.choice(range(256), size=3))) for i in range(len(phams))]
        pham_color_dict = dict(zip(phams,rgb_values))
        id = '|'.join([self.accessions,self.orf.id[0]])
        loc = np.where(pham_df.member == id)[0]
        color = pham_color_dict[list(pham_df.iloc[loc,0])[0]]
    def toJSON(self):
        return json.dumps(self, default=lambda o: o.__dict__,
            sort_keys=True, indent=4)

def cluster(phages, working_path):
    f = open(os.path.join(working_path, 'faa', 'orfs_pool.faa'), 'w')
    for phage in phages:
        for orf in phage.orfs:
            f.write('\n>{}|{}\n'.format(phage.accessions, orf.id[0]))
            f.write(orf.getSeq()[0])
    f.close()
    binb = '/usr/local/bin/'
    input_file = os.path.join(working_path, 'faa', 'orfs_pool.faa')
    createdb = ['{}/mmseqs'.format(binb),
                'createdb',
                input_file,
                os.path.join(working_path, 'cluster_out', 'DB')]
    print(createdb)
    subprocess.run(createdb, shell=False)
    cluster = ['{}/mmseqs'.format(binb),
               'cluster',
               '-v',
               '0',
               os.path.join(working_path, 'cluster_out', 'DB'),
               os.path.join(working_path, 'cluster_out', 'DB_clu'),
               os.path.join(working_path, 'tmp')]
    print(cluster)
    subprocess.run(cluster, shell=False)
    #mmseqs createtsv DB DB DB_clu DB_clu.tsv
    createtsv = ['{}/mmseqs'.format(binb),
                 'createtsv',
                 os.path.join(working_path,'cluster_out', 'DB'),
                 os.path.join(working_path,'cluster_out', 'DB'),
                 os.path.join(working_path,'cluster_out', 'DB_clu'),
                 os.path.join(working_path,'cluster_data','DB_clu.tsv')]
    print(createtsv)
    subprocess.run(createtsv, shell=False)
    return pd.read_csv(os.path.join(working_path, 'cluster_data/' 'DB_clu.tsv'), sep='\t',header=None, names=['representative','member'])

def blastn(phages, working_path):
    for phage in phages:
        output_name = os.path.join(working_path, 'fna', '{}.fna'.format(phage.accessions))
        SeqIO.convert(phage.filename, "genbank", output_name, "fasta")
    with open(os.path.join(working_path, 'commands.txt'), 'w') as handle:
        for query, target in list(combinations(glob.glob(os.path.join(working_path, 'fna', '*.fna')), 2)):
            q_name, t_name = os.path.basename(query).replace('.fna',''), os.path.basename(target).replace('.fna','')
            out = '{}_vs_{}.out'.format(q_name,t_name)
            blastx_cline = NcbiblastnCommandline(query=query, subject=target, outfmt=7,out=os.path.join(working_path, 'blast_out', out))
            handle.write(str(blastx_cline))
            handle.write('\n')
        handle.close()
    #parallel = ['/home/jbd_apps/bin/parallel', '<', os.path.join(working_path, 'commands.txt')]
    with open(os.path.join(working_path, 'commands.txt'), 'r') as f:
        print(f)
        subprocess.run(['/usr/local/bin/parallel'], stdin=f, check=True)
    results_di = {}
    for blast_out in glob.glob(os.path.join(working_path, "blast_out", '*.out')):
        results = pd.read_csv(blast_out, sep='\t',comment='#', names=['query', 'subject', 'identity', 'alignment' 'length', 'mismatches', 'gap_opens', 'q_start', 'q_end', 's_start', 's_end', 'evalue', 'bit_score'])
        results_di[os.path.basename(blast_out)] = results
    return results_di

def load_phages(phage_list):
    genome_list = []
    for phage in phage_list:
        with open(phage, "r") as input_handle:
            for record in SeqIO.parse(input_handle, "genbank"):
                locus = Locus(record.annotations, record.seq, phage)
                for i, feature in enumerate(record.features):
                    if feature.type == 'CDS' and 'protein_id' in feature.qualifiers.keys():
                        locus.add_feature(Orf(feature,i))
                    if feature.type == 'tRNA':
                        locus.add_feature(TRNA(feature))
                    if feature.type == 'repeat_region':
                        locus.add_feature(RepeatRegion(feature))
                genome_list.append(locus)
    return genome_list

def draw_shape(start,stop,strand, z, h, firstorf, lastorf):
    start, stop = int(start), int(stop)
    if strand == 1:
        x=(start, start+50, start, stop-50, stop, stop-50, start)
        y=(z, z+h/2, z+h, z+h, z+h/2, z, z)
    elif strand == -1:
        x=(start+50, start, start+50, stop, stop-50, stop, start+50)
        y=(z, z-h/2, z-h, z-h, z-h/2, z, z)
    else:
        z = z+0.5 #offset height
        x=(start, start, stop, stop, start)
        y=(z, z+h, z+h, z, z)
    return x,y

def draw_trna(start,stop,strand,z,h,firstorf,lastorf):
    start,stop = int(start), int(stop)
    if strand == 1:
        x=(start, start, stop, stop, start)
        y=(z, z+h, z+h, z, z)
    else:
        x=(start, start, stop, stop, start)
        y=(z, z-h, z-h, z, z)
    return x,y

def draw_repeat(start,stop,strand,z,h,firstorf,lastorf):
    start,stop = int(start), int(stop)
    if strand == 1:
        x=(start,start,stop,stop,start)
        y=(z,z+h,z+h,z,z)
    else:
        x=(start,start,stop,stop,start)
        y=(z,z-h,z-h,z,z)
    return x,y

def graphing(pham_df, phages, blast_di, order):
    phams = set(pham_df.representative)
    rgb_values = ['rgb{}'.format(tuple(np.random.choice(range(256), size=3))) for i in range(len(phams))]
    pham_color_dict = dict(zip(phams,rgb_values))
    #sorted_phages = sorted(phages, key=lambda x: x.annotations['organism'])
    order = [phages[x] for x in order]
    order_reduced = list(zip([x.accessions for x in order], range(len(order))))

    fig = go.Figure(layout={'width':1200,'height':1200})
    for z, phage in enumerate(order):
        fig.add_trace(go.Scatter(x=(0,len(phage.fna)),y=(z,z), mode='lines', line=dict(color='black', width=4, dash='dash')))
        shade = phage.load_syn_trace(blast_di, order_reduced, phages, z)
        [fig.add_trace(x) for x in shade]
    for z, phage in enumerate(order):
        try:
            [fig.add_trace(x) for x in phage.load_repeat_trace(z, 0.2)]
        except ValueError:
            pass
        if len(phage.orfs) > 0:
            [fig.add_trace(x) for x in phage.load_orf_trace(pham_color_dict, pham_df, z, 0.2)]
        if len(phage.tRNAs) > 0:
            [fig.add_trace(x) for x in phage.load_trna_trace(z, 0.2)]

    fig.update_layout(
        yaxis = dict(
            showgrid = False,
            zeroline = True,
        ),
        xaxis = dict(
            showgrid = True,
            zeroline = True,
            gridcolor = 'gray')

        )
    labels = [x.annotations['organism'] for x in order]
    fig.layout.plot_bgcolor = 'white'
    fig.layout.paper_bgcolor = 'white'
    fig.update_layout(showlegend=False)
    fig.update_layout(
        yaxis = dict(
            tickmode = 'array',
            tickvals = list(range(len(order))),
            ticktext = labels,
            )
        )

    fig.update_layout(
            height=(len(order)*50)+200,
            )
    fig.update_xaxes(range=[0, 8000])
    return fig

def defaultFig():
    fig = go.Figure()
    fig.layout.plot_bgcolor = 'white'
    fig.layout.paper_bgcolor = 'white'
    fig.update_layout(showlegend=False)
    fig.update_layout(
            yaxis = dict(
                visible = False),
            xaxis = dict(
                visible = False))
    return fig

def layout():
    layout = dbc.Container([
        html.H1('phamlite'),
        html.Div([
             dcc.Upload(
            id='upload-data',
            children=html.Div([
                'Drag and Drop or ',
                html.A('select genbank files')
            ]),
            style={
                'width': '100%',
                'height': '60px',
                'lineHeight': '60px',
                'borderWidth': '1px',
                'borderStyle': 'dashed',
                'borderRadius': '5px',
                'textAlign': 'center',
                'margin': '10px'
            },
            # Allow multiple files to be uploaded
            multiple=True
                ),
             html.Div(id='hidden_pham_df', style={'display': 'none'}),
             html.Div(id='hidden_phages', style={'display': 'none'}),
             html.Div(id='hidden_blast_di', style={'display': 'none'}),
             dcc.Dropdown(
                id = 'dropdown',
                options = [{'label':'select genomes','value': 0}],
                multi=True,
                ),
             dbc.Spinner(dcc.Graph(id='phamlite',
                figure = defaultFig(),
                 ),)
            ])
        ],className='mt-4')
    return layout
#if __name__ == '__main__':
app.layout = layout()


@app.callback([Output('dropdown','options'),
           Output('dropdown','value'),
           Output('hidden_pham_df','children'),
           Output('hidden_phages','children'),
           Output('hidden_blast_di','children')],
          [Input('upload-data','contents')],
          [State('upload-data','filename'),
           State('upload-data','last_modified')])
def load_dropdown(list_of_contents, list_of_names, list_of_dates):
    def randomString(stringLength=8):
        letters = string.ascii_lowercase
        return ''.join(random.choice(letters) for i in range(stringLength))
    storage_path = '/Users/matt/Desktop/phamlite_storage'
    if list_of_contents is not None:
        f  =  os.path.join(storage_path, randomString(8))
        print(f, flush=True)
        os.makedirs(f) #make working directory randomized string
        makedir(f) #make subdirectories in wd
        for i, contents in enumerate(list_of_contents):
            try:
                content_type, content_string = contents.split(',')
                file_content = base64.b64decode(content_string)
                with open( os.path.join(f,'{}.gb'.format(str(i)) ) ,"w+") as handle:
                    handle.write(file_content.decode("utf-8"))
            except ValueError:
                continue

        phage_list = glob.glob(os.path.join(f, '*.gb'))
        phages = load_phages(phage_list)
        print(phages, flush=True)
        blast_di = blastn(phages, f)
        print('CLUSTERING')
        pham_df = cluster(phages, f)
        print('done')
        labels = [ {'label':x.annotations['organism'], 'value':i} for i, x in enumerate(phages) ]
        pham_df = pham_df.to_json()
        phages = jsonpickle.encode(phages)
        blast_di = jsonpickle.encode(blast_di)
        return labels, list(range(len(labels))), pham_df, phages, blast_di
    else:
        print('empty', flush=True)
        raise dash.exceptions.PreventUpdate

@app.callback(Output('phamlite', 'figure'),
              [Input('dropdown', 'value'),
               Input('hidden_pham_df','children'),
               Input('hidden_phages','children'),
               Input('hidden_blast_di','children')])
def update_output(selected_order, pham_df, phages, blast_di):
    #pham_df = pd.read_json(pham_df)
    #phages = jsonpickle.decode(phages)
    #blast_di = jsonpickle.decode(blast_di)
    #fig = graphing(pham_df, phages, blast_di, selected_order)
    fig = {}
    return fig

app.run_server(host="0.0.0.0", port="8050")
