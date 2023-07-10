def contains_string_entire_column(df, substring):
    # cari string di seluruh kolom dataframe
    return df[df.apply(lambda row: row.astype(str).str.contains(substring, case=False).any(), axis=1)]

def contains_string_entire_column_boolean(df, substring):
    # cari string di seluruh kolom dataframe
    return df.apply(lambda row: row.astype(str).str.contains(substring, case=False).any(), axis=1)


def minmax(data):
    return (data - data.min())/ (data.max() - data.min())

def std_scale(data):
    return (data - data.mean()) / data.std()