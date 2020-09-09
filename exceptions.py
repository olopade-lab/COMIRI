#FIX THIS: override str function??

class IncorrectPathError(Exception):
    def __init__(self, path):
        self.message = ("Error: the following file path <{path}> "
                        "is incorrect.").format(path=path)

class WrongArgumentError(Exception):
    pass

class MissingArgumentError(Exception):
    def __init__(self, arg):
        self.message = ("Error: the following argument <{arg}> "
                        "is missing.").format(arg=arg)

class IncorrectInputFiles(Exception):
    pass