def tokenize(src: str):
    src = ''.join(src.split())  # remove whitespaces
    tokens = []
    reserve = ''

    def add_tok(tok: str):
        if tok != '':
            tokens.append(tok)

    def retrieve_store():
        nonlocal reserve
        tmp = reserve
        reserve = ''
        return tmp

    for c in src:
        if c.isalpha():
            if not reserve.isalpha() or c.isupper():
                add_tok(retrieve_store())
        elif c.isnumeric():
            if not reserve.isnumeric():
                add_tok(retrieve_store())
        else:
            add_tok(retrieve_store())
        reserve += c

    add_tok(retrieve_store())  # last token in reserve

    return tokens

