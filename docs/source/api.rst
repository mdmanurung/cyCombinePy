API reference
=============

This page is generated from the docstrings of every public symbol in
:mod:`cycombinepy`. Click through to each function for the full parameter
list, return type, and source link.

.. currentmodule:: cycombinepy

Preprocessing
-------------

.. autosummary::
   :toctree: generated/
   :nosignatures:

   transform_asinh
   normalize

Clustering
----------

.. autosummary::
   :toctree: generated/
   :nosignatures:

   create_som

Correction
----------

.. autosummary::
   :toctree: generated/
   :nosignatures:

   batch_correct
   correct_data
   CombatCorrectionError
   ConfoundedDesignError

Correction reports and constants
--------------------------------

.. currentmodule:: cycombinepy.correct

.. autosummary::
   :toctree: generated/
   :nosignatures:

   CORRECTED_LAYER
   CORRECTION_REPORT_UNS_KEY
   CORRECTION_REPORT_SCHEMA_VERSION

.. currentmodule:: cycombinepy

Evaluation
----------

.. autosummary::
   :toctree: generated/
   :nosignatures:

   compute_emd
   evaluate_emd
   compute_mad
   evaluate_mad

Detection
---------

.. autosummary::
   :toctree: generated/
   :nosignatures:

   detect_batch_effect
   detect_batch_effect_express

Utilities
---------

.. autosummary::
   :toctree: generated/
   :nosignatures:

   get_markers
   check_confound

I/O
---

.. currentmodule:: cycombinepy.io

.. autosummary::
   :toctree: generated/
   :nosignatures:

   read_fcs_dir

Plotting
--------

.. currentmodule:: cycombinepy.plotting

.. autosummary::
   :toctree: generated/
   :nosignatures:

   plot_density
   plot_dimred
   plot_emd_heatmap

Advanced evaluation
-------------------

.. currentmodule:: cycombinepy.evaluate

.. autosummary::
   :toctree: generated/
   :nosignatures:

   scib_metrics
