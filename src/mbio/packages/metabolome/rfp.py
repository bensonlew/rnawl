from sklearn.ensemble import RandomForestRegressor
import numpy as np
import pandas as pd
from sklearn.metrics import r2_score

def rf_pred(df, df_train, na_col):
    clf = RandomForestRegressor()
    fit = clf.fit(df_train.values.T, df.values)
    pred = clf.predict(df_train.values.T)
    for index,col in enumerate(df.index):
        if col in na_col:
            df.at[col] = pred[index]
    prob = r2_score(df.values, pred)
    return df,prob

def rf_permu(df, df_train, na_col, count_size=1):
    df,prob = rf_pred(df, df_train, na_col)
    # print "Predict score: %s" % prob
    if prob > 0.9 or count_size >= 10:
        print "Done in count_size %s, with r2 %s" % (count_size, prob)
        # print "Done in count_size %s" % count_size
        # print "check df:"
        # if type(df) == type(pd.Series()):
        #     print "check pass"
        # else:
        #     print "check error"
        return df
    else:
        df = rf_permu(df, df_train, na_col, count_size+1)
    return df

def get_value(df, method):
    df = df.replace(0, np.nan)
    if method=="median":
        value = df.median()
    return value

def missrf(df, df_train):
    if 0 in df.values:
        na_col = df[df==0].index
        value = get_value(df, "median")
        df = df.replace(0, value)
        index_list = df_train.index.tolist()
        index_list.remove(df.name)
        df_train = df_train[df_train.index.isin(index_list)]
        df = rf_permu(df, df_train, na_col)
    return df
'''
file = ""
data = pd.read(file, index_col=["metab_id"])
data.apply(missrf, axis=1, args=(data,))
'''