# -*- coding: utf-8 -*-
'''
This module defines :class:`RecordingChannelGroup`, a container for multiple
data channels.

:class:`RecordingChannelGroup` derives from :class:`BaseNeo`, from
:module:`neo.core.baseneo`.
'''

# needed for python 3 compatibility
from __future__ import absolute_import, division, print_function

import numpy as np

from neo.core.baseneo import BaseNeo


class RecordingChannelGroup(BaseNeo):
    '''
    A container for multiple data channels.

    This container have sereval purpose:
      * Grouping all :class:`AnalogSignalArray` inside a :class:`Block`
        across :class:`Segment`
      * Grouping :class:`RecordingChannel` inside a :class:`Block`. This
        case is *many to many* relation. It mean that a
        :class:`RecordingChannel` can belong to several group. A typical use
        case is tetrode (4 X :class:`RecordingChannel` inside a
        :class:`RecordingChannelGroup`).
      * Container of  :class:`Unit`. A neuron decharge (:class:`Unit`)
        can be seen by several electrodes (4 in tetrode case).

    *Usage 1* multi :class:`Segment` recording with 2 electrode array::

        >>> from neo.core import (Block, Segment, RecordingChannelGroup,
        ...                       AnalogSignalArray)
        >>> from quantities import nA, kHz
        >>> import numpy as np
        >>>
        >>> # create a Block with 3 Segment and 2 RecordingChannelGroup objects
        ,,, blk = Block()
        >>> for ind in range(3):
        ...     seg = Segment(name='segment %d' % ind, index=ind)
        ...     blk.segments.append(seg)
        ...
        >>> for ind in range(2):
        ...     rcg = RecordingChannelGroup(name='Array probe %d' % ind,
        ...                                 channel_indexes=np.arange(64))
        ...     blk.recordingchannelgroups.append(rcg)
        ...
        >>> # Populate the Block with AnalogSignalArray objects
        ... for seg in blk.segments:
        ...     for rcg in blk.recordingchannelgroups:
        ...         a = AnalogSignalArray(np.random.randn(10000, 64)*nA,
        ...                               sampling_rate=10*kHz)
        ...         rcg.analogsignalarrays.append(a)
        ...         seg.analogsignalarrays.append(a)

    *Usage 2* grouping channel::

        >>> from neo.core import (Block, RecordingChannelGroup,
        ...                       RecordingChannel)
        >>> import numpy as np
        >>>
        >>> # Create a Block
        ,,, blk = Block()
        >>>
        >>> # Create a new RecordingChannelGroup and add it to the Block
        ... rcg = RecordingChannelGroup(channel_names=np.array(['ch0',
        ...                                                     'ch1',
        ...                                                     'ch2']))
        >>> rcg.channel_indexes = np.array([0, 1, 2])
        >>> blk.recordingchannelgroups.append(rcg)
        >>>
        >>> # Create 3 RecordingChannel objects and add them to the Block
        ... for ind in range(3):
        ...     chan = RecordingChannel(index=ind)
        ...     rcg.recordingchannels.append(chan)  # <- many to many
        ,,,                                         # relationship
        ...     chan.recordingchannelgroups.append(rcg)  # <- many to many
        ...                                              #    relationship

    *Usage 3* dealing with :class:`Unit` objects::

        >>> from neo.core import Block, RecordingChannelGroup, Unit
        >>>
        >>> # Create a Block
        ... blk = Block()
        >>>
        >>> # Create a new RecordingChannelGroup and add it to the Block
        ... rcg = RecordingChannelGroup(name='octotrode A')
        >>> blk.recordingchannelgroups.append(rcg)
        >>>
        >>> # create several Unit objects and add them to the
        >>> # RecordingChannelGroup
        ... for ind in range(5):
        ...     unit = Unit(name = 'unit %d' % ind, description=
        ...                 'after a long and hard spike sorting')
        ...     rcg.units.append(unit)

    *Required attributes/properties*:
        None

    *Recommended attributes/properties*:
        :name: (str) A label for the dataset.
        :description: (str) Text description.
        :file_origin: (str) Filesystem path or URL of the original data file.
        :channel_names: (numpy.array 1D dtype='S')
            Names for each :class:`RecordingChannel`.
        :channel_indexes: (numpy.array 1D dtype='i')
            Index of each :class:`RecordingChannel`.

    Note: Any other additional arguments are assumed to be user-specific
            metadata and stored in :attr:`annotations`.

    *Container of*:
        :class:`RecordingChannel`
        :class:`AnalogSignalArray`
        :class:`Unit`

    '''

    def __init__(self, channel_names=None, channel_indexes=None, name=None,
                 description=None, file_origin=None, **annotations):
        '''
        Initialize a new :class:`RecordingChannelGroup` instance.
        '''
        # Inherited initialization
        # Sets universally recommended attributes, and places all others
        # in annotations
        BaseNeo.__init__(self, name=name, file_origin=file_origin,
                         description=description, **annotations)

        # Defaults
        if channel_indexes is None:
            channel_indexes = np.array([], dtype=np.int)
        if channel_names is None:
            channel_names = np.array([], dtype='S')

        # Store recommended attributes
        self.channel_names = channel_names
        self.channel_indexes = channel_indexes

        # Initialize containers for child objects
        self.analogsignalarrays = []
        self.units = []
        # Many to many relationship
        self.recordingchannels = []

        self.block = None

    def merge(self, other):
        '''
        Merge the contents of another RecordingChannelGroup into this one.

        For each :class:`RecordingChannel` in the other RecordingChannelGroup,
        if its name matches that of a :class:`RecordingChannel` in this block,
        the two RecordingChannels will be merged, otherwise it will be added as
        a new RecordingChannel. The equivalent procedure is then applied to
        each :class:`Unit`.

        For each array-type object in the other :class:`RecordingChannelGroup`,
        if its name matches that of an object of the same type in this segment,
        the two arrays will be joined by concatenation.
        '''
        for container in ("recordingchannels", "units"):
            lookup = dict((obj.name, obj) for obj in getattr(self, container))
            for obj in getattr(other, container):
                if obj.name in lookup:
                    lookup[obj.name].merge(obj)
                else:
                    getattr(self, container).append(obj)
        for container in ("analogsignalarrays",):
            objs = getattr(self, container)
            lookup = dict((obj.name, i) for i, obj in enumerate(objs))
            for obj in getattr(other, container):
                if obj.name in lookup:
                    ind = lookup[obj.name]
                    try:
                        newobj = getattr(self, container)[ind].merge(obj)
                    except AttributeError as e:
                        raise AttributeError("%s. container=%s, obj.name=%s, \
                                              shape=%s" % (e, container,
                                                           obj.name,
                                                           obj.shape))
                    getattr(self, container)[ind] = newobj
                else:
                    lookup[obj.name] = obj
                    getattr(self, container).append(obj)
        # TODO: merge annotations
