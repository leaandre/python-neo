# -*- coding: utf-8 -*-
"""

Class for reading data from Alpha Omega .map files.

This class is an experimental reader with important limitations.
See the source code for details of the limitations.
The code of this reader is of alpha quality and received very limited testing.

This code is written from the incomplete file specifications available in:

[1] AlphaMap Data Acquisition System User's Manual Version 10.1.1
Section 5 APPENDIX B: ALPHAMAP FILE STRUCTURE, pages 120-140
Edited by ALPHA OMEGA Home Office: P.O. Box 810, Nazareth Illit 17105, Israel
http://www.alphaomega-eng.com/

and from the source code of a C software for conversion of .map files to
.eeg elan software files :

[2] alphamap2eeg 1.0, 12/03/03, Anne CHEYLUS - CNRS ISC UMR 5015

Supported : Read

@author : sgarcia, Florent Jaillet

"""

# NOTE: For some specific types of comments, the following convention is used:
# "TODO:" Desirable future evolution
# "WARNING:" Information about code that is based on broken or missing
# specifications and that might be wrong


# Main limitations of this reader:
# - The reader is only able to load data stored in data blocks of type 5
#   (data block for one channel). In particular it means that it doesn't
#   support signals stored in blocks of type 7 (data block for multiple
#   channels).
#   For more details on these data blocks types, see 5.4.1 and 5.4.2 p 127 in
#   [1].
# - Rather than supporting all the neo objects types that could be extracted
#   from the file, all read data are returned in AnalogSignal objects, even for
#   digital channels or channels containing spiking informations.
# - Digital channels are not converted to events or events array as they
#   should.
# - Loading multichannel signals as AnalogSignalArrays is not supported.
# - Many data or metadata that are avalaible in the file and that could be
#   represented in some way in the neo model are not extracted. In particular
#   scaling of the data and extraction of the units of the signals are not
#   supported.
# - It received very limited testing, exlusively using python 2.6.6. In
#   particular it has not been tested using Python 3.x.
#
# These limitations are mainly due to the following reasons:
# - Incomplete, unclear and in some places innacurate specifications of the
#   format in [1].
# - Lack of test files containing all the types of data blocks of interest
#   (in particular no file with type 7 data block for multiple channels where
#   available when writing this code).
# - Lack of knowledge of the Alphamap software and the associated data models.
# - Lack of time (especially as the specifications are incomplete, a lot of
#   reverse engineering and testing is required, which makes the development
#   of this IO very painful and long).


# needed for python 3 compatibility
from __future__ import (absolute_import, division)

# specific imports
import datetime
import logging
import os
import struct


# file no longer exists in Python3
try:
    file
except NameError:
    import io
    file = io.BufferedReader

# note neo.core need only numpy and quantities
import numpy as np
import quantities as pq

from neo.io.baseio import BaseIO
from neo.core import Block, Segment, AnalogSignal, EventArray, SpikeTrain
from neo.io.tools import create_many_to_one_relationship
from neo.io.tools import populate_RecordingChannel


class AlphaOmegaIO(BaseIO):
    """
    Class for reading data from Alpha Omega .map files (experimental)

    This class is an experimental reader with important limitations.
    See the source code for details of the limitations.
    The code of this reader is of alpha quality and received very limited
    testing.

    Usage:
        >>> from neo import io
        >>> r = io.AlphaOmegaIO(filename = 'File_AlphaOmega_1.map')
        >>> blck = r.read_block(lazy = False, cascade = True)
        >>> print blck.segments[0].analogsignals
    """

#    Attributes:
#        - option_readLevel (string):
#            2 option :
#                - 'AnalogSignal' :
#                   DEFAULT. level signal are loaded in AnalogSignal neo object
#                - 'SpikeTrain' :
#                   level signal are loaded in SpikeTrain neo object
#        - option_readWaveform (boolean):
#            if True, wavaforms are read for SpikeTrain object
#
#    Methods:
#        Public:
#            - option_readLevel
#            - option_readWaveform
#            - read_block
#        Private:
#            - _read_analogData
#            - _read_digitalData
#            - _read_spikeData
#            - _read_portData
#            - _annotate_block
#            - _annotate_block2
#            - _annotate_blockb
#            - _count_samples
#            - _select_availableChannel
#
#
#    Inherits from:
#        neo.io.BaseIO

    is_readable = True  # This is a reading only class
    is_writable = False  # writting is not supported

    # This class is able to directly or inderectly read the following kind of
    # objects
    supported_objects = [Block, Segment, AnalogSignal, EventArray, SpikeTrain]
    # TODO: Add support for other objects that should be extractable from .map
    # files (AnalogSignalArray, Event, Epoch?, Epoch Array?, Spike?)

    readable_objects = [Block]  # This class can only return a Block
    # TODO: create readers for different type of objects (Segment,
    # AnalogSignal,...)

    writeable_objects = []  # This class is not able to write objects

    # This is for GUI stuff : a definition for parameters when reading.
    read_params = {Block: []}

    write_params = None  # Writing is not supported, so no GUI stuff

    name = 'AlphaOmega'
    extensions = ['map', 'mpx']
    mode = 'file'

    # write is not supported so I do not overload write method from BaseIO

    def __init__(self, filename=None):
        """
        Arguments:
            filename (string): The Alpha Omega file name. DEFAULT = None.

        Other Attributes:
            option_readLevel (string): name of neo object to store level datas
            option_readWaveform (boolean): if True, waveforms are read for
                SpikeTrain object

        """
        BaseIO.__init__(self)
        self.filename = filename
        self.option_readLevel = "AnalogSignal"
        self.option_readWaveform = False

        # Logger definition
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.WARNING)
        ch = logging.StreamHandler()
        ch.setLevel(logging.WARNING)
        formatter = logging.Formatter('%(name)s - %(levelname)s - %(message)s')
        ch.setFormatter(formatter)
        self.logger.addHandler(ch)

#==============================================================================
#   Public options
#==============================================================================

    def option_readLevel(self, option_readLevel):
        """
        option_readLevel setter.

        Arguments:
            option_readLevel (string): name of neo object to store level datas
                options :
                        - 'AnalogSignal' : DEFAULT.
                        - 'SpikeTrain'

        """
        options = ['AnalogSignal', 'SpikeTrain']

        try:
            option_readLevel in options
        except:
            raise Exception("Option %s is not available." % (option_readLevel))
        else:
            self.option_readLevel = option_readLevel
#==============================================================================
#   Private read data method
#==============================================================================

    def _read_analogData(self, fid, list_blocks, file_blocks):
        """
        Read datas present in block 5, and return datas and informations

        Data type read:
            - Continuous Data
            - Level Data
        
        Arguments:
            fid (file): file object, ready to read
            list_blocks (list): list of indexes of data blocks
            file_blocks (list): list of blocks

        Returns:
            signal_array (np.array): array fill with the analog signal
            start_signal (int): start of analog signal in absolute time

        """
        # NOTE: we assume that the data block are in time range, so the first
        # block contains the start_time and with the last block, we can find
        # the stop_time for signal.

        # Start time
        count = self._count_samples(file_blocks[list_blocks[0]]['m_length'])
        fid.seek(file_blocks[list_blocks[0]]['pos'] + 6 + (count * 2))
        start_signal = struct.unpack('I', fid.read(4))[0]

        # End time
        count = self._count_samples(
            file_blocks[list_blocks[len(list_blocks) - 1]]['m_length'])
        fid.seek(file_blocks[list_blocks[len(list_blocks) - 1]]['pos'] +
                 6 + (count * 2))
        end_signal = struct.unpack('I', fid.read(4))[0] + count

        # Initialize signal array
        chan_signal_len = end_signal - start_signal
        signal_array = np.empty(chan_signal_len)
        signal_array.fill(np.nan)

        for ind_block in list_blocks:
            # Count number of samples
            count = self._count_samples(file_blocks[ind_block]['m_length'])

            # Find start of block
            fid.seek(file_blocks[ind_block]['pos'] + 6 + (count * 2))
            start = struct.unpack('I', fid.read(4))[0]

            # Data reading
            start_block = start - start_signal
            end_block = start - start_signal + count
            fid.seek(file_blocks[ind_block]['pos'] + 6)
            signal_array[start_block:end_block] = np.fromfile(
                fid,
                dtype=np.int16,
                count=count)

        return signal_array, start_signal

    def _read_digitalData(self, fid, list_blocks, file_blocks, chan_len,
                          ind_chan, ind):
        """
        Read datas present in block 5, and return 2 temp_array

        Data type read:
            - Digital Data

        Arguments:
            fid (file): file object, ready to read
            list_blocks (list): list of indexes of data blocks
            file_blocks (list): list of blocks
            chan_len (np.array): array of len(total data) for a channel
            ind_chan (int): index of channel in 'chan_len' vector
            ind (int): index in the data vector

        Returns:
            temp_array_up (np.array): array fill with time for UP event
            temp_array_down (np.array): array fill with time for DOWN event

        """
        # Initialize array
        temp_array = np.empty(chan_len[ind_chan], dtype=np.float64)
        labels = np.empty(chan_len[ind_chan], dtype=np.int16)

        for ind_block in list_blocks:
            # Count number of samples
            count = self._count_samples(file_blocks[ind_block]['m_length'])

            # Position in file for times datas
            fid.seek(file_blocks[ind_block]['pos'] + 6)
            temp_array[ind:ind + count] = np.fromfile(
                fid,
                dtype=np.int32,
                count=count)

            # Position in file for labels datas
            fid.seek(file_blocks[ind_block]['pos'] + 6 + 4)
            labels[ind:ind + count] = np.fromfile(
                fid,
                dtype=np.int16,
                count=count)
            ind += count

        temp_array *= pq.ms

        # Separate in Up and Down temp_array
        temp_array_up = temp_array[labels == 1]
        temp_array_down = temp_array[labels == 0]

        return temp_array_up, temp_array_down

    def _read_spikeData(self, fid, list_blocks, file_blocks):
        """
        Read datas present in block 5, and return objects for SpikeTrain
        creation

        Data type read:
            - Level Data

        Arguments:
            fid (file): file object, ready to read
            list_blocks (list): list of indexes of data blocks
            file_blocks (list): list of blocks

        Returns:
            times (np.array): array fill with time (float)
            t_start (int) : start time of data
            t_stop (int) : stop time of data
            waveforms (Quantity 3D) : object containing spike waveform

        """
        waveforms = None

        # Start time
        count = self._count_samples(file_blocks[list_blocks[0]]['m_length'])
        fid.seek(file_blocks[list_blocks[0]]['pos'] + 6 + (count * 2))
        t_start = struct.unpack('I', fid.read(4))[0]

        # End time
        count = self._count_samples(
            file_blocks[list_blocks[len(list_blocks) - 1]]['m_length'])
        fid.seek(file_blocks[list_blocks[len(list_blocks) - 1]]['pos'] +
                 6 + (count * 2))
        t_stop = struct.unpack('I', fid.read(4))[0] + count

        # Initialize times and waveforms array
        count = self._count_samples(file_blocks[list_blocks[0]]['m_length'])
        times = np.empty(len(list_blocks))
        times.fill(np.nan)
        if self.option_readWaveform:
            waveforms = np.empty((len(list_blocks), 1, count))

        for ind, ind_block in enumerate(list_blocks):
            # Count number of samples
            count = self._count_samples(file_blocks[ind_block]['m_length'])

            # Times data reading
            fid.seek(file_blocks[ind_block]['pos'] + 6 + (count * 2))
            time = struct.unpack('I', fid.read(4))[0]
            times[ind] = time

            # Waveforms data reading
            if self.option_readWaveform:
                # Initialize spike array
                spike_array = np.empty(count)
                spike_array.fill(np.nan)
                # Waveform reading
                fid.seek(file_blocks[ind_block]['pos'] + 6)
                spike_array[0:count] = np.fromfile(
                    fid,
                    dtype=np.int16,
                    count=count)
                waveforms[ind,:,:].fill(np.nan)
                waveforms[ind,:,:] = spike_array

        return times, t_start, t_stop, waveforms

    def _read_portData(self, fid, list_blocks, file_blocks, chan_len,
                       ind_chan, ind):
        """
        Read datas present in block 5, and return 2 temp_array

        Data type read:
            - Port Data

        Arguments:
            fid (file): file object, ready to read
            list_blocks (list): list of indexes of data blocks
            file_blocks (list): list of blocks
            chan_len (np.array): array of len(total data) for a channel
            ind_chan (int): index of channel in 'chan_len' vector
            ind (int): index in the data vector

        Returns:
            temp_array (np.array): array fill with time (float)
            labels (np.array): array fill with labels (int)

        """
        # Initialize time and label array
        # NOTE: array are normally crate with np.empty method, but it cause
        # a problem for neomatlab reading. It actually work with np.zeros
        # method.
        temp_array = np.zeros(chan_len[ind_chan], dtype=np.float64)
        labels = np.zeros(chan_len[ind_chan], dtype=np.int16)

        for ind_block in list_blocks:
            # Count number of samples
            count = self._count_samples(file_blocks[ind_block]['m_length'])

            # Position in file for events or labels datas
            fid.seek(file_blocks[ind_block]['pos'] + 6)
            labels[ind:ind + count] = np.fromfile(
                fid,
                dtype=np.int16,
                count=count)

            # Position in file for times datas
            fid.seek(file_blocks[ind_block]['pos'] + 6 + 2)
            temp_array[ind:ind + count] = np.fromfile(
                fid,
                dtype=np.int32,
                count=count)
            ind += count

        return temp_array, labels

#==============================================================================
#   Private annotate block method
#==============================================================================

    def _annotate_block(self, neo_object, block):
        """
        Annotate neo_object with data present in information block

        Block supported:
            - 2
            - b

        Arguments:
            neo_object (neo_object): EventArray, AnalogSignal, ...
            block (dict): dictionnary which contains block informations

        Returns:
            neo_object (neo_object): neo_object with annotations

        """
        block_available = ['2','b']

        try:
            block['m_TypeBlock'] in block_available
        except:
            raise Exception("Annotate is not available for block type %s"
                % (block['m_TypeBlock']))
        else:
            if block['m_TypeBlock'] == '2':
                self._annotate_block2(neo_object, block)
            elif block['m_TypeBlock'] == 'b':
                self._annotate_blockb(neo_object, block)

        return neo_object

    def _annotate_block2(self, neo_object, block):
        """
        Annotate neo_object with data present in block 2

        Arguments:
            neo_object (neo_object): EventArray, AnalogSignal, ...
            block (dict): dictionnary which contains block informations

        Returns:
            neo_object (neo_object): neo_object with annotations

        """
        neo_object.annotate(isInput=block['m_isInput'])
        neo_object.annotate(numColor=block['m_numColor'])
        neo_object.annotate(channel_type=block['type_subblock'])

        if block['type_subblock'] != 'digital':
            neo_object.annotate(Amplitude=block['m_Amplitude'] * pq.volt)
            
        if block['type_subblock'] == 'digital':
            neo_object.annotate(SaveTrigger=block['m_SaveTrigger'])
            neo_object.annotate(Duration=block['m_Duration'] * pq.ms)
            neo_object.annotate(PreviousStatus=block['m_PreviousStatus'])

        if block['type_subblock'] == 'level':
            neo_object.annotate(LevelValue=block['m_LevelValue'] * pq.volt)
            neo_object.annotate(TrgMode=block['m_TrgMode'])
            neo_object.annotate(YesRMS=block['m_YesRMS'])
            neo_object.annotate(bAutoScale=block['m_bAutoScale'])

        if (block['type_subblock'] == 'level' or
            block['type_subblock'] == 'external_trigger'):
                neo_object.annotate(nSpikeCount=block['m_nSpikeCount'])
                neo_object.annotate(
                    nPreTrigmSec=block['m_nPreTrigmSec'] * pq.ms)
                neo_object.annotate(
                    nPostTrigmSec=block['m_nPostTrigmSec'] * pq.ms)

        if (block['type_subblock'] == 'external_trigger' or
            block['type_subblock'] == 'digital'):
                neo_object.annotate(channel_index=int(block['m_numChannel']))

        return neo_object

    def _annotate_blockb(self, neo_object, block):
        """
        Annotate neo_object with data present in block b

        Arguments:
            neo_object (neo_object): EventArray, AnalogSignal, ...
            block (dict): dictionnary which contains block informations

        Returns:
            neo_object (neo_object): neo_object with annotations

        """
        neo_object.annotate(BoardNumber=block['m_BoardNumber'])
        neo_object.annotate(channel_index=block['m_numChannel'])
        neo_object.annotate(channel_type='port')

        return neo_object

#==============================================================================
#   Private specific method
#==============================================================================

    def _count_samples(self, m_length):
        """
        Count the number of signal samples available in a type 5 data block

        Arguments:
            m_length (int): length of block 5

        Returns:
            count (int): number of samples in block 5

        """
        # INFOS: for information about type 5 data block, see [1].
        # -6 corresponds to the header of block 5, and the -2 take into
        # account the fact that last 2 values are not available as the 4
        # corresponding bytes are coding the time stamp of the beginning
        # of the block
        count = int(((m_length - 6) / 2) - 2)
        return count

    def _select_availableChannel(self, groups, list_available_channel):
        """
        Sort on groups and subgroups with available channel.
        See read_block() Step 1b for groups construction structure.

        Arguments:
            groups (list): list of group (dict)
            list_available_channel (list) : list of available channel

        Returns:
            groups (list): list of available group (dict)

        """

        # Step 1: delete unavailable channels
        for group in groups:
            for subgroup in group['subgroups']:
                for channel in subgroup['channels']:
                    if channel not in list_available_channel:
                        subgroup['channels'].remove(channel)

        # Step 2: delete empty subgroups (without available channels)
        for group in groups:
            for subgroup in group['subgroups']:
                if len(subgroup['channels']) == 0:
                    group['subgroups'].remove(subgroup)

        # Step 3: delete empty groups (without available subgroups)
        for group in groups:
            if len(group['subgroups']) == 0:
                groups.remove(group)

        return groups

#==============================================================================
#   Public read method
#==============================================================================

    def read_block(self, lazy=False, cascade=True):
                   # the 2 first keyword arguments are imposed by neo.io API
        """
        Return a Block.

        Arguments:
            lazy (bool):
                if True, the numerical data is not loaded, only
                properties and metadata
            cascade (bool):
                if False, only a single object is loaded, without its
                references to other objects

        Returns:
            blck (Block)

        """

        # create the neo Block that will be returned at the end
        blck = Block(file_origin=os.path.basename(self.filename))
        blck.file_origin = os.path.basename(self.filename)

        fid = open(self.filename, 'rb')

        #=====================================================================
        # Step 1: read the headers of all the data blocks to load the file
        # structure
        #=====================================================================
        # NOTE: in the following, the word "block" is used in the sense
        # used in the alpha-omega specifications (ie a data chunk in the
        # file), rather than in the sense of the usual Block object in neo

        pos_block = 0  # position of the current block in the file
        file_blocks = []  # list of data blocks available in the file
        list_type_block = []  # list of data blocks type in the file

        informations = {}  # dict containing general informations
        groups = []  # list containing groups informations
        out_infos = [  # infos not returned in informations dict
            "pos",
            "blank",
            "blank2",
            "m_length",
            "m_TypeBlock",
            "m_nextBlock",
            "m_EraseCount",
            "m_Reserved",
            "m_placeMainWindow"]
        DAP_infos = [  # infos for log
            "DAP_buffers_filling",
            "RMS_value",
            "num_channel",
            "sample_count"]
        i_group = 0  # index of group
        i_subgroup = 0  # index of subgroup

        if not cascade:  # we read only the main header

            m_length, m_TypeBlock = struct.unpack('Hcx', fid.read(4))
            # m_TypeBlock should be 'h', as we read the first block
            blck = HeaderReader(
                fid,
                dict_header_type.get(m_TypeBlock, Type_Unknown)).read_f()
            blck.update({
                'm_length': m_length,
                'm_TypeBlock': m_TypeBlock,
                'pos': pos_block})
            file_blocks.append(blck)
            list_type_block.append(m_TypeBlock)

            # Recup date and time informations
            if m_TypeBlock == 'h':
                blck.rec_datetime = datetime.datetime(
                    blck['m_date_year'],
                    blck['m_date_month'],
                    blck['m_date_day'],
                    blck['m_time_hour'],
                    blck['m_time_minute'],
                    blck['m_time_second'],
                    10000 * blck['m_time_hsecond'])
                    # the 10000 is here to convert m_time_hsecond from
                    # centisecond to microsecond

        else:  # cascade == True

            seg = Segment(file_origin=os.path.basename(self.filename))
            seg.file_origin = os.path.basename(self.filename)
            blck.segments.append(seg)

            while True:
                first_4_bytes = fid.read(4)

                if len(first_4_bytes) < 4:
                    # we have reached the end of the file
                    break
                else:
                    m_length, m_TypeBlock = struct.unpack('Hcx', first_4_bytes)

                block = HeaderReader(
                    fid,
                    dict_header_type.get(m_TypeBlock, Type_Unknown)).read_f()
                block.update({
                    'm_length': m_length,
                    'm_TypeBlock': m_TypeBlock,
                    'pos': pos_block})
                list_type_block.append(m_TypeBlock)

                #========================================
                # a - Read of subblock for block 2 and 7
                #========================================

                if m_TypeBlock == '2':
                    # The beggining of the block of type '2' is identical for
                    # all types of channels, but the following part depends on
                    # the type of channel. So we need a special case here.

                    # WARNING: How to check the type of channel is not
                    # described in the documentation. So here I use what is
                    # proposed in the C code [2].
                    # According to this C code, it seems that the 'm_isAnalog'
                    # is used to distinguished analog and digital channels, 
                    # and 'm_Mode' encodes the type of analog channel:
                    # 0 for continuous, 1 for level, 2 for external trigger.
                    # But in some files, I found channels that seemed to be
                    # continuous channels with 'm_Modes' = 128 or 192. So I
                    # decided to consider every channel with 'm_Modes'
                    # different from 1 or 2 as continuous. I also couldn't
                    # check that values of 1 and 2 are really for level and
                    # external trigger as I had no test files containing data
                    # of this types.

                    type_subblock = 'unknown_channel_type'
                    description = Type2_SubBlockUnknownChannels
                    block.update({'m_Name': 'unknown_name'})

                    if block['m_isAnalog'] == 0:
                        # digital channel
                        type_subblock = 'digital'
                        description = Type2_SubBlockDigitalChannels
                    elif block['m_isAnalog'] == 1:
                        # analog channel
                        # NOTE : analog channel block have a different
                        # header than digital channel. Here is the reading of
                        # the identical part for analog channel
                        description = Type2_SubBlockAnalogChannels
                        subblock = HeaderReader(fid, description).read_f()
                        block.update(subblock)

                        if block['m_Mode'] == 1:
                            # level channel
                            type_subblock = 'level'
                            description = Type2_SubBlockLevelChannels
                        elif block['m_Mode'] == 2:
                            # external trigger channel
                            type_subblock = 'external_trigger'
                            description = Type2_SubBlockExtTriggerChannels
                        else:
                            # continuous channel
                            type_subblock = 'continuous (Mode %i)' % (
                                block['m_Mode'])
                            description = Type2_SubBlockContinuousChannels

                    subblock = HeaderReader(fid, description).read_f()
                    block.update(subblock)
                    block.update({'type_subblock': type_subblock})

                if m_TypeBlock == '7':
                    # The beggining of the block '7' is identical, datas
                    # next are different according to the documentation [1].

                    description = Type7_DataSubblockUnknown

                    if block['FINT'] == -111:
                        description = Type7_DataSubblockDAPBuffers
                    elif block['FINT'] == -222:
                        description = Type7_DataSubblockRMS
                    elif block['FINT'] == -444:
                        description = Type7_DataSubblockRestart
                    elif block['FINT'] == -333:
                        description = Type7_DataSubblockDataLoss
                    elif block['FINT'] == -557 or block['FINT'] == -558:
                        description = Type7_DataSubblockStopAnalogOutput

                    subblock = HeaderReader(fid, description).read_f()
                    block.update(subblock)

                #=============================================================
                #  b - Annotate and log informations
                #=============================================================
                #  INFOS: groups informations structure:
                #       groups = [group, group, group, ...]
                #       group = {
                #           'name': group_name,
                #           'number': group_number,
                #           'Z_Order': appearence_order,
                #           'subgroups': [subgroup, subgroups, subgroup, ...]}
                #       subgroup = {
                #           'name': group_name,
                #           'number': group_number,
                #           'Z_Order': appearence_order,
                #           'TypeOverLap': typeoverlap,
                #           'channels': [channel, channel, channel, ...]}
                #       general and groups informations are stored in
                #       annotations of Block.

                # Recup date and time informations
                elif m_TypeBlock == 'h':
                    blck.rec_datetime = datetime.datetime(
                        block['m_date_year'],
                        block['m_date_month'],
                        block['m_date_day'],
                        block['m_time_hour'],
                        block['m_time_minute'],
                        block['m_time_second'],
                        10000 * block['m_time_hsecond'])
                        # the 10000 is here to convert m_time_hsecond from
                        # centisecond to microsecond

                    seg.rec_datetime = blck.rec_datetime.replace()
                    # I couldn't find a simple copy function for datetime,
                    # using replace without arguments is a twisted way to make
                    # a copy

                # Recup general informations
                if m_TypeBlock in ['h', '0']:
                    for infos in block:
                        if (infos in out_infos or
                            infos.startswith("m_MainWindow_") or
                            infos.startswith("m_date_") or
                            infos.startswith("m_time_")):
                                pass
                        else:
                            if infos in Type1_quantities:
                                informations[infos[2:]] = (
                                    block[infos] * Type1_quantities[infos])
                            else:
                                informations[infos[2:]] = block[infos]
                                
                if block['m_TypeBlock'] == '1':
                    for infos in block:
                        if infos not in out_infos:
                            if infos in Type1_quantities:
                                informations["Board_%i_%s" % (
                                 block["m_Number"],
                                 infos.strip('m_'))] = (
                                     block[infos] * Type1_quantities[infos])
                            else:
                                informations["Board_%i_%s" % (
                                 block["m_Number"],
                                 infos.strip('m_'))] = block[infos]

                # Recup groups informations
                if m_TypeBlock == '3':
                    dict_group = {
                        'name': block['m_nameGroup'],
                        'number': block['m_Number'],
                        'Z_Order': block['m_Z_Order'],
                        'subgroups': []}
                    groups.append(dict_group)
                    i_group += 1
                    i_subgroup = 0

                if m_TypeBlock == '4':
                    i_subgroup += 1
                    dict_subgroup = {
                        'name': block['m_Name'],
                        'number': block['m_Number'],
                        'Z_Order': block['m_Z_Order'],
                        'TypeOverlap': block['m_TypeOverlap'],
                        'channels': []}
                    groups[i_group - 1]['subgroups'].append(dict_subgroup)
                    for info in block:
                        if (info.startswith('m_numChannel') and
                            block[info] != 0):
                                groups[i_group - 1]['subgroups'] \
                                      [i_subgroup - 1]['channels'] \
                                      .append(block[info])

                # Raise errors about DAP informations
                if m_TypeBlock == '7':
                    error = "\nBoard %i - Error %i" % (
                        block['m_numBoard'],
                        block['FINT'])
                    for info in block:
                        if info in DAP_infos:
                            error += " - %s = %s" % (info, block[info])
                    if block['FINT'] == -111:
                        error = "DAP Buffers - " + error
                    elif block['FINT'] == -222:
                        error = "RMS - " + error
                    elif block['FINT'] == -444:
                        error = "DAP Restart - " + error
                    elif block['FINT'] == -333:
                        error = "Data Loss " + error
                        error += " - Loss between %i - %i" % (
                            block["first_lost_sample"],
                            block["last_lost_sample"])
                    elif block['FINT'] == -666:
                        error = "Start Analog Output - " + error
                    elif block['FINT'] == -557 or block['FINT'] == -558:
                        error = "Stop Analog Output - " + error

                    self.logger.info(error)

                file_blocks.append(block)
                pos_block += m_length
                fid.seek(pos_block)

            #==============================================================
            # Step 2: find the available channels
            #==============================================================
            #       NOTE: Block 2 contains informations for continuous,
            #       digital, level and trigger channels. Block b contains
            #       informations for Port channels, except the sampling_rate
            #       present in corresponding digital channels (number 11033 
            #       to 110040)

            list_chan = []  # list containing indexes of channel blocks
            ind_port = None  # indice for port information contain in block 2

            for ind_block, block in enumerate(file_blocks):
                if block['m_TypeBlock'] == '2' or block['m_TypeBlock'] == 'b':
                    list_chan.append(ind_block)
                    if block['m_numChannel'] == 11033:
                        ind_port = ind_block
                    

            #================================================================
            # Step 3: find blocks containing data for the available channels
            #================================================================

            list_data = []  # list of lists of indexes of data blocks
                            # corresponding to each channel
            list_available_channel = []  # list of available numChannel

            for ind_chan, chan in enumerate(list_chan):
                list_data.append([])
                num_chan = file_blocks[chan]['m_numChannel']
                for ind_block, block in enumerate(file_blocks):
                    if block['m_TypeBlock'] == '5':
                        list_available_channel.append(block['m_numChannel'])
                        if block['m_numChannel'] == num_chan:
                            list_data[ind_chan].append(ind_block)

            #================================================================
            # Step 4: compute the length (number of samples) of the channels
            #================================================================

            chan_len = np.zeros(len(list_data), dtype=np.int)

            for ind_chan, list_blocks in enumerate(list_data):
                for ind_block in list_blocks:
                    chan_len[ind_chan] += self._count_samples(
                        file_blocks[ind_block]['m_length'])

            #====================================================
            # Step 5: find channels for which data are available
            #====================================================

            ind_valid_chan = np.nonzero(chan_len)[0]

            #=================================================================
            # Step 6: load the data
            #     TODO give the possibility to load data as AnalogSignalArrays
            #=================================================================

            for ind_chan in ind_valid_chan:
                ind = 0  # index in the data vector
                list_blocks = list_data[ind_chan]

                # read time stamp for the beginning of the signal
                form = '<l'  # reading format
                ind_block = list_blocks[0]

                count = self._count_samples(file_blocks[ind_block]['m_length'])
                fid.seek(file_blocks[ind_block]['pos'] + 6 + (count * 2))
                buf = fid.read(struct.calcsize(form))
                val = struct.unpack(form, buf)
                start_index = val[0]

                # find the signal type
                numChannel = file_blocks[ind_block]['m_numChannel']
                type_signal = ''
                for block in file_blocks:
                    if (block['m_TypeBlock'] == '2' and
                        block['m_numChannel'] == numChannel):
                            type_signal = block['type_subblock']
                    elif (block['m_TypeBlock'] == 'b' and
                          block['m_numChannel'] == numChannel):
                            type_signal = 'port'

                sampling_rate = \
                    file_blocks[list_chan[ind_chan]]['m_SampleRate']
                if int(sampling_rate) != 0:
                    t_start = start_index / sampling_rate
                else:
                    t_start = start_index

                #=============================================================
                # CONTINUOUS signal and LEVEL signal with ANALOGSIGNAL option
                #=============================================================
                if ('continuous' in type_signal or
                    (type_signal == 'level' and
                     self.option_readLevel == 'AnalogSignal')):

                    if lazy:
                        ana_sig = AnalogSignal(
                            np.array([]) * pq.dimensionless,
                            sampling_rate=sampling_rate * pq.kHz,
                            t_start=t_start * pq.ms,
                            name=file_blocks[list_chan[ind_chan]]['m_Name'],
                            file_origin=os.path.basename(self.filename),
                            channel_index = int(file_blocks \
                                [list_chan[ind_chan]]['m_numChannel']),
                            units=pq.dimensionless)
                        ana_sig.lazy_shape = chan_len[ind_chan]

                    else:
                        (signal_array, t_start) = self._read_analogData(
                            fid,
                            list_blocks,
                            file_blocks)

                        t_start = (t_start / sampling_rate)
                        amplitude = file_blocks[
                            list_chan[ind_chan]]['m_Amplitude']
                        quantum = float(amplitude / signal_array.max())
                        signal_array *= quantum

                        ana_sig = AnalogSignal(
                            signal_array * pq.volt,
                            sampling_rate=sampling_rate * pq.kHz,
                            t_start=t_start * pq.ms,
                            name=file_blocks[list_chan[ind_chan]]['m_Name'],
                            file_origin=os.path.basename(self.filename),
                            channel_index = int(file_blocks \
                                [list_chan[ind_chan]]['m_numChannel']),
                            units=pq.volt)

                    ana_sig = self._annotate_block(
                        ana_sig,
                        file_blocks[list_chan[ind_chan]])

                    seg.analogsignals.append(ana_sig)

                #=====================================
                # LEVEL signal with SPIKETRAIN option
                #=====================================
                if (type_signal == 'level' and
                     self.option_readLevel == 'SpikeTrain'):

                    if lazy:
                        spike_train = SpikeTrain(
                            np.array([]) * pq.ms,
                            sampling_rate=sampling_rate * pq.kHz,
                            t_start=t_start * pq.ms,
                            t_stop=None,
                            left_sweep=file_blocks[list_chan[ind_chan]] \
                                ['m_nPreTrigmSec'] * pq.ms,
                            name=file_blocks[list_chan[ind_chan]]['m_Name'],
                            file_origin=os.path.basename(self.filename),
                            units=pq.ms)
                        spike_train.lazy_shape = chan_len[ind_chan]

                    else:
                        (times,
                         t_start,
                         t_stop,
                         waveforms) = self._read_spikeData(
                            fid,
                            list_blocks,
                            file_blocks)

                        times /= sampling_rate
                        t_start = (t_start / sampling_rate)
                        t_stop = (t_stop / sampling_rate)

                        spike_train = SpikeTrain(
                            times * pq.ms,
                            sampling_rate=sampling_rate * pq.kHz,
                            t_start=t_start * pq.ms,
                            t_stop=t_stop * pq.ms,
                            left_sweep=file_blocks[list_chan[ind_chan]] \
                                ['m_nPreTrigmSec'] * pq.ms,
                            name=file_blocks[list_chan[ind_chan]]['m_Name'],
                            file_origin= os.path.basename(self.filename),
                            units=pq.ms)
                        if self.option_readWaveform:
#                            amplitude = file_blocks[
#                            list_chan[ind_chan]]['m_Amplitude']
#                            quantum = amplitude / (2**((
#                               signal_array.min() * -1) + signal_array.max()) -1)
#                            signal_array *= quantum
                            spike_train.waveforms = waveforms * pq.dimensionless

                    spike_train.channel_index = int(file_blocks \
                        [list_chan[ind_chan]]['m_numChannel'])
                    spike_train = self._annotate_block(
                        spike_train,
                        file_blocks[list_chan[ind_chan]])
                    seg.spiketrains.append(spike_train)

                #==================
                #  DIGITAL signal
                #==================
                elif type_signal == 'digital':
                    if lazy:
                        dig_sig_up = EventArray(
                            np.array([]) * pq.dimensionless,
                            labels=np.array([], dtype='S'),
                            name=file_blocks[list_chan[ind_chan]]['m_Name'] \
                                + "_Up",
                            file_origin=os.path.basename(self.filename),
                            units=pq.dimensionless)
                        dig_sig_down = EventArray(
                            np.array([]) * pq.dimensionless,
                            labels=np.array([], dtype='S'),
                            name=file_blocks[list_chan[ind_chan]]['m_Name'] \
                                + "_Down",
                            file_origin=os.path.basename(self.filename))

                    else:
                        (temp_array_up,
                         temp_array_down) = self._read_digitalData(
                            fid,
                            list_blocks,
                            file_blocks,
                            chan_len,
                            ind_chan,
                            ind)
                        # Treatment of temp_array
                        temp_array_up /= sampling_rate
                        temp_array_down /= sampling_rate

                        labels_up = np.array(temp_array_up.shape[0],
                                             dtype='S')
                        labels_down = np.array(temp_array_down.shape[0],
                                               dtype='S')

                        dig_sig_up = EventArray(
                            temp_array_up * pq.ms,
                            labels=labels_up,
                            name=file_blocks[list_chan[ind_chan]]['m_Name'] \
                                + "_Up",
                            file_origin = os.path.basename(self.filename))
                        dig_sig_down = EventArray(
                            temp_array_down * pq.ms,
                            labels=labels_down,
                            name=file_blocks[list_chan[ind_chan]]['m_Name'] \
                                + "_Down",
                            file_origin = os.path.basename(self.filename))

                    dig_sig_up = self._annotate_block(
                        dig_sig_up,
                        file_blocks[list_chan[ind_chan]])
                    dig_sig_down = self._annotate_block(
                        dig_sig_down,
                        file_blocks[list_chan[ind_chan]])

                    dig_sig_up.annotate(sampling_rate=sampling_rate * pq.kHz)
                    dig_sig_down.annotate(sampling_rate=sampling_rate * pq.kHz)

                    seg.eventarrays.append(dig_sig_up)
                    seg.eventarrays.append(dig_sig_down)


                #==================
                #  PORT signal
                #==================
                elif type_signal == 'port':

                    if lazy:
                        port_sig = EventArray(
                            np.array([]) * pq.dimensionless,
                            labels=np.array([], dtype='S'),
                            name=file_blocks[list_chan[ind_chan]]['m_Name'],
                            file_origin=os.path.basename(self.filename))

                    else:
                        temp_array, labels = self._read_portData(
                            fid,
                            list_blocks,
                            file_blocks,
                            chan_len,
                            ind_chan,
                            ind)
                        # WARNING: sampling_rate is not present in
                        # the block b which refer to the channel (value of 0.0)
                        # the sampling_rate is present in block 2 for 
                        # corresponding digital channels
                        sampling_rate = file_blocks[
                            list_chan[ind_port]]['m_SampleRate']
                        temp_array /= sampling_rate # Treatment of temp_array
                        t_start /= sampling_rate
                        
                        port_sig = EventArray(
                            temp_array * pq.ms,
                            labels=labels,
                            name=file_blocks[list_chan[ind_chan]]['m_Name'],
                            file_origin = os.path.basename(self.filename))

                    port_sig = self._annotate_block(
                        port_sig,
                        file_blocks[list_chan[ind_chan]])
                    port_sig.annotate(sampling_rate=sampling_rate * pq.kHz)
                    port_sig.annotate(t_start=t_start * pq.ms)
                    seg.eventarrays.append(port_sig)

        fid.close()

        #==============================================
        # Step 7: recup informations and annotate them
        #==============================================

        # Select group and subgroup with available data
        groups = self._select_availableChannel(groups,list_available_channel)

        # Annotations
        for info in informations:
            blck.annotations[info] = informations[info]
        blck.annotate(groups=groups)

        # Specific annotations for tools (To put out after)
#        blck.annotate(list_block=list_type_block)
#        blck.annotate(file_blocks=file_blocks)

        if cascade:
            populate_RecordingChannel(blck, remove_from_annotation=False)
#            create_many_to_one_relationship(blck)

        return blck




"""
Information for special types in [1]:

_dostime_t type definition:
struct dos_time_t
{
 unsigned char hour; /* hours (0-23)*/
 unsigned char minute; /* minutes (0-59)*/
 unsigned char second; /* seconds (0-59) */
 unsigned char hsecond; /* seconds/ 100 (0-99)*/
}

_dosdate_t type definition:
struct _dosdate_t
{
 unsigned char day;       /* day of month( 1-31) */
 unsigned char month;     /* month (1-12) */
 unsigned int year;       /* year (1980-2099) */
 unsigned char dayofweek; /* day of week (0 = Sunday) */
}

WINDOWPLACEMENT16 type definition (according to WINE source code):
typedef struct
{
    UINT16   length;
    UINT16   flags;
    UINT16   showCmd;
    POINT16  ptMinPosition;
    POINT16  ptMaxPosition;
    RECT16   rcNormalPosition;
} WINDOWPLACEMENT16,*LPNONCLIENTMETRICS16;

POINT16 struct
{
    INT16 x;
    INT16 y;
}

RECT16 Struct
{
    INT16 bottom;
    INT16 left;
    INT16 right;
    INT16 top;
}

"""

max_string_len = '80s'

# maximal length of variable length strings in the file
# WARNING: I don't know what is the real value here. According to [1] p 139
# it seems that it could be 20. Some tests would be needed to check this.

# WARNING: A cleaner way to handle strings reading is suitable. Currently I
# read a buffer of max_string_len bytes and look for the C "end of string"
# character ('\x00'). It would be better either to read characters until
# reaching '\x00' or to read the exact number of characters needed, if the
# length of a string can be deduced from the lentgh of the block and the number
# of bytes already read (it seems possible, at least for certain block types).

# The name of the keys in the folowing dicts are chosen to match as closely as
# possible the names in document [1]

TypeH_Header = [
    ('m_nextBlock','l'),
    ('m_version','h'),
    ('m_time_hour','B'),
    ('m_time_minute','B'),
    ('m_time_second','B'),
    ('m_time_hsecond','B'),
    ('m_date_day','B'),
    ('m_date_month','B'),
    ('m_date_year','H'),
    ('m_date_dayofweek','B'),
    ('blank', 'x'), # one byte blank for 2 bytes alignement
    ('m_minTime','d'),
    ('m_maxTime','d'),
    ('m_EraseCount','l'),
    ('m_mapVersion','b'),
    ('m_ApplicationName','10s'),
    ('m_ResourceVersion','4s'), # WARNING: present only in version 65,
    ('blank2','x'),
    ('m_Reserved','l')] # WARNING: Mpx version. Must be checked


Type0_SetBoards = [
    ('m_nextBlock','l'),
    ('m_BoardCount','h'),
    ('m_GroupCount','h'),
    ('m_MainWindow_length','H'),
    ('m_MainWindow_flags','H'),
    ('m_MainWindow_showCmd','H'),
    ('m_MainWindow_ptMinPosition_x','H'),
    ('m_MainWindow_ptMinPosition_y','H'),
    ('m_MainWindow_ptMaxPosition_x','H'),
    ('m_MainWindow_ptMaxPosition_y','H'),
    ('m_MainWindow_rcNormalPosition_bottom','H'),
    ('m_MainWindow_rcNormalPosition_left','H'),
    ('m_MainWindow_rcNormalPosition_right','H'),
    ('m_MainWindow_rcNormalPosition_top','H')] # WARNING: order unknown

Type1_Boards = [ # WARNING: needs to be checked
    ('m_nextBlock','l'),
    ('m_Number','h'),
    ('m_countChannel','h'),
    ('m_countAnIn','h'),
    ('m_countAnOut','h'),
    ('m_countDigIn','h'),
    ('m_countDigOut','h'),
    ('m_TrigCount', 'h'), # not defined in 5.3.3 but appears in 5.5.1 and
                          # seems to really exist in files
    # WARNING: check why 'm_TrigCount is not in the C code [2]
    ('m_Amplitude','f'),
    ('m_cSampleRate','f'), # sample rate seems to be given in kHz
    ('m_Duration','f'),
    ('m_nPreTrigmSec','f'),
    ('m_nPostTrigmSec','f'),
    ('m_TrgMode','h'),
    ('m_LevelValue','h'), # after this line, 5.3.3 is wrong,
                          # check example in 5.5.1 for the right fields
    # WARNING: check why the following part is not corrected in the C code [2]
    ('m_nSamples','h'),
    ('m_LevelFactorRMS','f'),
    ('m_ScaleFactor','f'),
    ('m_DapTime','f'),
    ('m_nameBoard', '35s'),
    ('m_DiscMaxValue','b'), # WARNING: should this exist?
    ('m_DiscMinValue','b'), # WARNING: should this exist?
    ('blank','x')] # WARNING : may be before m_DiscMaxValue or after.

Type2_DefBlocksChannels = [
    # common parameters for all types of channels
    ('m_nextBlock','l'),
    ('m_isAnalog','h'),
    ('m_isInput','h'),
    ('m_numChannel','h'),
    ('m_numColor','h')]

Type2_SubBlockAnalogChannels = [
    ('blank', '2x'), # WARNING : this is not in the specs but it seems needed
    ('m_Mode','h')]

Type2_SubBlockContinuousChannels = [
    # continuous channels parameters
    ('m_Amplitude','f'),
    ('m_SampleRate','f'),
    ('m_ContBlkSize','h'),
    ('m_ModeSpike','h'), # WARNING: the C code [2] uses usigned short here
    ('m_Duration','f'),
    ('m_bAutoScale','h'),
    ('m_Name', max_string_len)]

Type2_SubBlockLevelChannels = [
    # level channels parameters
    ('m_Amplitude','f'),
    ('m_SampleRate','f'),
    ('m_nSpikeCount','h'),
    ('m_ModeSpike','h'),
    ('m_nPreTrigmSec','f'),
    ('m_nPostTrigmSec','f'),
    ('m_LevelValue','h'),
    ('m_TrgMode','h'),
    ('m_YesRMS','h'),
    ('m_bAutoScale','h'),
    ('m_Name', max_string_len)]

Type2_SubBlockExtTriggerChannels = [ # WARNING: untested
    # external trigger channels parameters
    ('m_Amplitude','f'),
    ('m_SampleRate','f'),
    ('m_nSpikeCount','h'),
    ('m_ModeSpike','h'),
    ('m_nPreTrigmSec','f'),
    ('m_nPostTrigmSec','f'),
    ('m_TriggerNumber','h'),
    ('m_Name', max_string_len)]

Type2_SubBlockDigitalChannels = [
    # digital channels parameters
    ('m_Mode','h'),
    ('m_SampleRate','f'),
    ('m_SaveTrigger','h'),
    ('m_Duration','f'),
    ('m_PreviousStatus','h'), # WARNING: check difference with C code here
    ('m_Name', max_string_len)]

Type2_SubBlockUnknownChannels = [
    # WARNING: We have a mode that doesn't appear in our spec, so we don't
    # know what are the fields.
    # It seems that for non-digital channels the beginning is
    # similar to continuous channels. Let's hope we're right...
    ('blank', '2x'),
    ('m_Amplitude','f'),
    ('m_SampleRate','f')]
    # there are probably other fields after...

Type6_DefBlockTrigger = [
    ('m_nextBlock','l'),
    ('m_Number','h'),
    ('m_countChannel','h'),
    ('m_StateChannels','h'),
    ('m_numChannel1','h'),
    ('m_numChannel2','h'),
    ('m_numChannel3','h'),
    ('m_numChannel4','h'),
    ('m_numChannel5','h'),
    ('m_numChannel6','h'),
    ('m_numChannel7','h'),
    ('m_numChannel8','h'),
    ('m_Name',max_string_len)]

Type3_DefBlockGroup = [
    ('m_nextBlock','l'),
    ('m_Number','h'),
    ('m_Z_Order','h'),
    ('m_countSubGroups','h'),
    ('m_MainWindow_length','H'),
    ('m_MainWindow_flags','H'),
    ('m_MainWindow_showCmd','H'),
    ('m_MainWindow_ptMinPosition_x','H'),
    ('m_MainWindow_ptMinPosition_y','H'),
    ('m_MainWindow_ptMaxPosition_x','H'),
    ('m_MainWindow_ptMaxPosition_y','H'),
    ('m_MainWindow_rcNormalPosition_bottom','H'),
    ('m_MainWindow_rcNormalPosition_left','H'),
    ('m_MainWindow_rcNormalPosition_right','H'),
    ('m_MainWindow_rcNormalPosition_top','H'), # WARNING: order unknown
    ('m_NetLoc','h'),
    ('m_locatMax','4I'),
    ('m_nameGroup', max_string_len)] # 'c' in documentation

Type4_DefBlockSubgroup = [
    ('m_nextBlock','l'),
    ('m_Number','h'),
    ('m_TypeOverlap','h'),
    ('m_Z_Order','h'),
    ('m_countChannel','h'),
    ('m_NetLoc','h'),
    ('m_location','4I'),
    ('m_bIsMaximized','h'),
    ('m_numChannel1','H'),
    ('m_numChannel2','H'),
    ('m_numChannel3','H'),
    ('m_numChannel4','H'),
    ('m_numChannel5','H'),
    ('m_numChannel6','H'),
    ('m_numChannel7','H'),
    ('m_numChannel8','H'),
    ('m_Name', max_string_len)]

Type5_DataBlockOneChannel = [
    ('m_numChannel','h')]
    # WARNING: 'm_numChannel' (called 'm_Number' in 5.4.1 of [1]) is supposed
    # to be uint according to 5.4.1 but it seems to be a short in the files
    # (or should it be ushort ?)

# WARNING: In 5.1.1 page 121 of [1], they say "Note: 5 is used for demo
# purposes, 7 is used for real data", but looking at some real datafiles,
# it seems that block of type 5 are also used for real data...

Type7_DataBlockMultipleChannels = [
    ('m_numBoard', 'H'), # WARNING: unknown true type
    ('FINT','h')]

Type7_DataSubblockDAPBuffers = [
    # FINT = - 111
    ('DAP_buffers_filling','h')]

Type7_DataSubblockRMS = [ # WARNING: untested
    # FINT = - 222
    ('RMS_value', 'b')]

Type7_DataSubblockRestart = [ # WARNING: untested
    # FINT = - 444
    ('num_channel', 'h')]

Type7_DataSubblockDataLoss = [
    # FINT = - 333
    ('num_channel', 'H'),
    ('first_lost_sample','i'),
    ('last_lost_sample','i')]

Type7_DataSubblockStartAnalogOutput = [ # WARNING: untested
    # FINT = - 666
    ]

Type7_DataSubblockStopAnalogOutput = [ # WARNING: untested
    # FINT = - 557 or - 558
    ('sample_count','i')]

Type7_DataSubblockUnknown = [ # WARNING: untested
    ('value','h')]

TypeP_DefBlockPeriStimHist = [ # WARNING: present in version 100
    ('m_numChannel','h'),
    ('m_Position','4I'),
    ('m_isStatVisible','h'),
    ('m_DurationSec','f'),
    ('m_Rows','i'),
    ('m_DurationSecPre','f'),
    ('m_Bins','i'),
    ('m_NoTrigger','h')]

TypeF_DefBlockFRTachogram = [ # WARNING: untested
    ('m_numChannel','h'),
    ('m_Position','4I'),
    ('m_isStatVisible','h'),
    ('m_DurationSec','f'),
    ('m_AutoManualScale','i'),
    ('m_Max','i')]

TypeR_DefBlockRaster = [
    ('m_numChannel','h'),
    ('m_Position','4I'),
    ('m_isStatVisible','h'),
    ('m_DurationSec','f'),
    ('blank','2x'), # WARNING: position not sure
    ('m_Rows','i'),
    ('m_NoTrigger','h'),
    ('nextblock','Hcx')]

TypeI_DefBlockISIHist = [ # WARNING: untested
    ('m_numChannel','h'),
    ('m_Position','4I'),
    ('m_isStatVisible','h'),
    ('m_DurationSec','f'),
    ('m_Bins','i'),
    ('m_TypeScale','i')]

Type8_MarkerBlock = [ # WARNING: untested
    ('m_numChannel','h'),
    ('m_Time','l')] # WARNING: check what's the right type here.
    # It seems that the size of time_t type depends on the system typedef,
    # I put long here but I couldn't check if it is the right type

Type9_ScaleBlock = [ # WARNING: untested
    ('m_numChannel','h'),
    ('m_Scale','f')]

Typeb_DigPortDef = [
    ('m_BoardNumber','i'),
    ('m_numChannel','h'),
    ('m_SampleRate','f'),
#    ('m_PrevValue','H'),    # WARNING : seems absent in file
    ('m_Name', max_string_len)]

TypeS_DefStreamData = [ # WARNING: untested
    ('m_nextBlock','l'),
    ('m_Number','h'),
    ('m_SampleRate','f'),
    ('m_Name', max_string_len)]

TypeE_StreamDataBlock = [ # WARNING: untested
    ('m_TypeEvent','c'),
    ('m_uTimeStamp','L'),
    ('m_StreamData', max_string_len)]

Type_Unknown = []


dict_header_type = {
                    'h' : TypeH_Header,
                    '0' : Type0_SetBoards,
                    '1' : Type1_Boards,
                    '2' : Type2_DefBlocksChannels,
                    '6' : Type6_DefBlockTrigger,
                    '3' : Type3_DefBlockGroup,
                    '4' : Type4_DefBlockSubgroup,
                    '5' : Type5_DataBlockOneChannel,
                    '7' : Type7_DataBlockMultipleChannels,
                    'P' : TypeP_DefBlockPeriStimHist,
                    'F' : TypeF_DefBlockFRTachogram,
                    'R' : TypeR_DefBlockRaster,
                    'I' : TypeI_DefBlockISIHist,
                    '8' : Type8_MarkerBlock,
                    '9' : Type9_ScaleBlock,
                    'b' : Typeb_DigPortDef,
                    'S' : TypeS_DefStreamData,
                    'E' : TypeE_StreamDataBlock,
                    }

Type1_quantities = {
    "m_minTime": pq.s,
    "m_maxTime": pq.s,
    "m_Amplitude": pq.volt,
    "m_cSampleRate": pq.kHz,
    "m_Duration": pq.ms,
    "m_nPreTrigmSec": pq.ms,
    "m_nPostTrigmSec": pq.ms,
    "m_LevelValue": pq.volt,
    }


class HeaderReader():
    '''
    Class for metadatas reading in AlphaOmega.map files

    '''

    def __init__(self, fid, description):
        '''
        Arguments:
            fid (file object) : file ready for reading
            description (dict) : dictionnary with :
                key (string) = value_name
                value (string) = C type format reading

        '''
        self.fid = fid
        self.description = description

    def read_f(self, offset=None):
        '''
        Method for reading file with unpack method

        Arguments:
            offset(int) : DEFAULT = None.

        Returns:
            d (dict): dictionnary with metadatas

        '''
        if offset is not None:
            self.fid.seek(offset)

        d = {}  # dict with metadatas
        for key, fmt in self.description:
            fmt = '<' + fmt # insures use of standard sizes
            buf = self.fid.read(struct.calcsize(fmt))
            # value reading
            if len(buf) != struct.calcsize(fmt) :
                return None
            val = list(struct.unpack(fmt, buf))
            # value treatment
            for i, ival in enumerate(val):
                if hasattr(ival, 'split'):
                    val[i] = ival.split('\x00', 1)[0]
            if len(val) == 1:
                val = val[0]
            # value storage
            d[key] = val
        return d
