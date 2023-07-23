#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from tsfresh import extract_features
df = pd.read_table('/Users/aikosekita/Documents/AD_longitudinal/longitudinal_set_timeseries_data.txt')
X = extract_features(df, column_id="patient_id", column_sort="date")
X.to_csv('/Users/aikosekita/Documents/AD_longitudinal/longitudinal_set_timeseries_extract_feature.csv', index=False)