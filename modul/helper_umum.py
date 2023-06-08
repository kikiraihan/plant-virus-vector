def contains_string_entire_column(df, substring):
    # cari string di seluruh kolom dataframe
    return df[df.apply(lambda row: row.astype(str).str.contains(substring, case=False).any(), axis=1)]

def contains_string_entire_column_boolean(df, substring):
    # cari string di seluruh kolom dataframe
    return df.apply(lambda row: row.astype(str).str.contains(substring, case=False).any(), axis=1)