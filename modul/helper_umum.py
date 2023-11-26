import json

def contains_string_entire_column(df, substring):
    # cari string di seluruh kolom dataframe
    return df[df.apply(lambda row: row.astype(str).str.contains(substring, case=False).any(), axis=1)]

def contains_set_of_strings_entire_column(df, substrings):
    # cari set string di seluruh kolom dataframe
    return df[df.apply(lambda row: any(substring.lower() in str(cell).lower() for substring in substrings for cell in row), axis=1)]

def contains_string_entire_column_boolean(df, substring):
    # cari string di seluruh kolom dataframe
    return df.apply(lambda row: row.astype(str).str.contains(substring, case=False).any(), axis=1)


def minmax(data):
    return (data - data.min())/ (data.max() - data.min())

def std_scale(data):
    return (data - data.mean()) / data.std()

def report_back(progress,message):
    data=json.dumps({
        "progress":progress,
        'message':message,
        })
    return f'data: {data}\n\n'

def pemecah_generator(ini):
    def main_generator():
        yield "Awal"
        result = yield from ini
        # yield result  # Menggunakan hasil yang dikembalikan dari sub-generator
        yield result # terakhir ambil return

    return [i for i in main_generator()][-1]