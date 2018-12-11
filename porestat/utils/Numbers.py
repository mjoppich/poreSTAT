def toFloat(sText, nDefault=None):
    """

    :param sText: a text string
    :param nDefault: a default value to return if sText is not a number
    :return: int value of sText or None if not a number
    """
    if type(sText) is int:
        return float(sText)

    if type(sText) is float:
        return float(sText)

    if (len(sText) == 0):
        return nDefault

    if str(sText[0]).isdigit():

        try:
            iValue = float(sText)
            return iValue

        except:

            return nDefault

    return nDefault


def toInt(sText, nDefault=None):
    """

    :param sText: a text string
    :param nDefault: a default value to return if sText is not a number
    :return: int value of sText or None if not a number
    """
    if type(sText) is int:
        return sText

    if type(sText) is float:
        return int(sText)

    if (len(sText) == 0):
        return nDefault

    if str(sText[0]).isdigit():

        try:
            iValue = int(sText)
            return iValue

        except:

            return nDefault

    return nDefault


def isNumber(sText):
    """

    :param sText: a text string
    :return: returns TRUE if sText is a number
    """
    if len(sText) == 0:
        return False

    if str(sText[0]).isdigit():

        try:
            iValue = int(sText)
            return True

        except:

            try:
                fValue = float(sText)
                return True
            except:

                return False

    return False


def toNumber(sText, nDefault=None):
    """

    :param sText: a text string
    :return: returns TRUE if sText is a number
    """
    oValue = toInt(sText, None)

    if oValue is None:
        return toFloat(sText, nDefault)

    return oValue
