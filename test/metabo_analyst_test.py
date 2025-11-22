import os
import pickle

import pandas as pd
import pytest

from src.metabolome.metabo_analyst import MetaboAnalyst


def test_get_metabolites():
    res = MetaboAnalyst.query(['Tryptophan', 'l-kynurenine'])
    print(res)




