class IncorrectPathError(Exception):
    def __init__(self, path):
        Exception.__init__(self, ("Error: the following file path <{path}> "
                                    "is incorrect.").format(path=path))

class WrongArgumentError(Exception):
    pass

class MissingArgumentError(Exception):
    def __init__(self, arg):
        Exception.__init__(self, ("Error: the following argument <{arg}> "
                                    "is missing.").format(arg=arg))
#override the str method

class SubprocessError(Exception):
    def __init__(self, program):
        Exception.__init__(self, ("An error occurred while running <{program}>, "
                                "the pipeline did not run to completion.").format(program=program))

class IncorrectInputFiles(Exception):
    pass