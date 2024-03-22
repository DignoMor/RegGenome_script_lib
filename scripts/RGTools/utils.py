
def str2bool(bool_str):
    '''
    Convert String to boolean value
    '''

    if not bool_str:
        return False

    if bool_str.upper() == "FALSE":
        return False

    if bool_str.upper() == "None":
        return False

    return True
