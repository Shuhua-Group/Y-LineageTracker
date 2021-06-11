import configparser
import os


def getConfig():

    config = configparser.ConfigParser()
    config.read(os.path.split(os.path.realpath(__file__))[0]+'/Config.ini')

    return config
