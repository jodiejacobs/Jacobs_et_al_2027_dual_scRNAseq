#!/usr/bin/env python3
"""Run DeepSig under Keras 3 by aliasing `keras` to the legacy tf_keras API."""
import sys, types, tf_keras
m = types.ModuleType("keras"); m.__path__ = []
for sub in ("models", "backend", "utils", "preprocessing"):
    if hasattr(tf_keras, sub):
        setattr(m, sub, getattr(tf_keras, sub)); sys.modules["keras." + sub] = getattr(tf_keras, sub)
sys.modules["keras"] = m
seq = types.ModuleType("keras.preprocessing.sequence")
seq.pad_sequences = tf_keras.utils.pad_sequences
sys.modules["keras.preprocessing.sequence"] = seq
from deepsig.deepsig import main
sys.exit(main())
