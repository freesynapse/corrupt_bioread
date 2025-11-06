#!/usr/bin/python3

#%%
import io
import math
import sys

from langsmith import expect
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import bioread_corrupt
from bioread_corrupt.reader import ChunkBuffer

# argument parsing
if len(sys.argv) != 3:
   print('Usage: read_corrupt_acq.py <ACQ FILE> <EXPORT CSV FILE>')
   sys.exit(1)
#
fname = sys.argv[1]
if fname[-4:] != '.acq':
   print('ERROR: invalid file extension of AcqKnowledge data file. Exit.')
   sys.exit(2)
# 
export_fname = sys.argv[2]
if export_fname[-4:] != '.csv':
   export_fname += '.csv'

# VALID = False

# if VALID:
    # fname = './data/valid.acq'
# else:
# fname = './data/corrupt.acq'

c = bioread_corrupt.read(fname)

# Even if event markers are truncated, the data is still in the file.
# Here, we read the data at the data start offset
# In the files I checked this is at byte 45590, might be different depending on the 
# version of AckNowledge
#
data_offset = c.channel_dtype_headers[0].offset


d = {}
var_names = [c_.name for c_ in c.channel_headers]
channel_idxs = list(range(len(var_names)))

print('Channel data:')
block_sz_bytes = 0
for i, var in zip(channel_idxs, var_names):
    d[i] = {}
    d[i]['descr'] = var
    d[i]['np_dtype'] = c.channel_dtype_headers[i].numpy_dtype
    d[i]['sz_bytes'] = 2 if d[i]['np_dtype'] == '<i2' else 8
    if d[i]['np_dtype'] == '<i2':
        block_sz_bytes += 2
    elif d[i]['np_dtype'] == '<f8':
        block_sz_bytes += 8
    else:
        raise TypeError
    
    d[i]['offset_bytes'] = c.channel_dtype_headers[i].offset
    d[i]['raw_data_scale'] = c.channels[i].raw_scale_factor
    d[i]['raw_data_offset'] = c.channels[i].raw_offset
    
    print(f'\t{var:32}: offset = {c.channel_dtype_headers[i].offset:5d}, sz = {d[i]['sz_bytes']}')

assert(block_sz_bytes == np.sum([c_.sample_size for c_ in c.channels]))
print(f'block_sz_bytes = {block_sz_bytes}')
# hack; the data starts 4 bytes after the last offset
data_offset = d[channel_idxs[-1]]['offset_bytes'] + 4

# Open file for reading and calculate the readable data size (i.e. the number of complete
# data time steps). End of file is corrupt, so assumes the rest of the file is data. 
#
file = open(fname, 'rb')
file_sz = file.seek(0, io.SEEK_END)
# data_sz = file_sz - d[0]['offset']
data_sz = file_sz - data_offset
block_count = math.floor(data_sz / block_sz_bytes)
# block_count = 784010
file.close()

# read data into numpy array
np_dtypes = np.dtype([(d[i]['descr'], d[i]['np_dtype']) for i in range(len(d.keys()))])

# read the whole data block at once (no chunking)
file = open(fname, 'rb')
# file.seek(d[0]['offset'])
file.seek(data_offset)
buf = file.read(block_count * block_sz_bytes)
data = np.frombuffer(buf, 'b')

# byte pattern per data block, corresponding to the index of the channel, 0 for SAP etc
data_block_byte_pattern = []
for i, var in enumerate(var_names):
    for _ in range(d[i]['sz_bytes']):
        data_block_byte_pattern.append(i)
data_block_byte_pattern = np.array(data_block_byte_pattern)

# create a mask with the respective data buffer, corresponding to the 
channel_idx_byte_mask = np.tile(data_block_byte_pattern, block_count)

# create buffers, one per channel, to hold the data
cbuffers = [ChunkBuffer(idx) for idx in channel_idxs]

# read and reinterpret raw bytes into buffers using the mask
for i in channel_idxs:
    buffer = cbuffers[i]
    buffer.raw_buffer = data[channel_idx_byte_mask == i]
    buffer.raw_buffer.dtype = d[i]['np_dtype']
    buffer.buffer = buffer.raw_buffer * d[i]['raw_data_scale'] + d[i]['raw_data_offset']

# valid data, compare to BioRead
#if VALID:
#    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(10, 4))
#    ax0.plot(c.channels[0].data[:200])
#    ax0.set_title('BioRead')
#    ax0.spines[['top', 'right']].set_visible(False)
#    ax1.plot(cbuffers[0].buffer[:200])
#    ax1.set_title('myown')
#    ax1.spines[['top', 'right']].set_visible(False)
#    plt.show()
#else:
#    fig, ax = plt.subplots(1)
#    ax.plot(cbuffers[0].buffer[:200])
#    ax.set_title('myown -- corrupt file!')
#    ax.spines[['top', 'right']].set_visible(False)
#    plt.show()
    
#%%

# csv export
#
print(f'Exporting to \'{export_fname}\'...')
cnames = var_names.copy()
cnames.insert(0, 'time_s')
cnames.insert(0, 'time_idx')
#
n = cbuffers[0].buffer.shape[0]
df = pd.DataFrame(data=None, columns=cnames)

for i, var in zip(channel_idxs, var_names):
    df[var] = cbuffers[i].buffer
df['time_idx'] = np.arange(0, df.shape[0])
df['time_s'] = np.arange(0, df.shape[0]) / float(c.channels[0].samples_per_second)

df.to_csv(export_fname, sep=',')

# sys.exit(0)

#%%
# checking the csv; sanity check
print('Sanity check: SAP plot.')
df1 = pd.read_csv('export.csv', sep=',')
fig, ax = plt.subplots(figsize=(12, 5))
ax.plot(df1['time_s'], df1['SAP'], lw=0.5, c='firebrick')
ax.set_xlim(0, df1['time_s'].to_numpy()[-1])
ax.set_title(f'SAP read from \'{fname}\'')
ax.set_xlabel('Time (s)')
ax.set_ylabel('SAP (mmHg)')
ax.spines[['top', 'right']].set_visible(False)
plt.show()
