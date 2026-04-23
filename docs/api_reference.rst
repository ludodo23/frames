API Reference
=============

FrameGraph
----------

.. doxygenclass:: frames::FrameGraph
   :project: frames
   :members:

Frame Strategies — Rotation
-----------------------------

.. doxygenstruct:: frames::BConstantRotation
   :project: frames
   :members:

.. doxygenclass:: frames::BFixedAtEpochRotation
   :project: frames
   :members:

.. doxygenstruct:: frames::BSampledRotation
   :project: frames
   :members:

Frame Strategies — Translation
--------------------------------

.. doxygenstruct:: frames::BConstantTranslation
   :project: frames
   :members:

.. doxygenclass:: frames::BFixedAtEpochTranslation
   :project: frames
   :members:

.. doxygenstruct:: frames::BSampledTranslation
   :project: frames
   :members:

Sampled Data
------------

.. doxygenstruct:: frames::SampledData
   :project: frames
   :members:

Eigen Backend
-------------

.. doxygenstruct:: frames::EigenBackend
   :project: frames
   :members:

.. doxygenenum:: frames::Axis
   :project: frames

.. doxygenstruct:: frames::Intrinsic
   :project: frames

.. doxygenstruct:: frames::Extrinsic
   :project: frames

.. doxygenfunction:: frames::eulerToQuaternion
   :project: frames
