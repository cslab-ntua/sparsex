#!/usr/bin/env python
# vim:expandtab:softtabstop=4:shiftwidth=4:tabstop=8

import traceback

def debug(func):
    """Decorator function. Prints the arguments and the return value of a function."""

    def quote(x):
        if type(x) == str: return '"%s"' % x[:100]
        else: return repr(x)[:300]

    def wrapper(*args, **kwargs):
        #traceback.print_tb()
        params = [quote(x) for x in args] + ["%s=%s" % (k, quote(v)) for k, v in kwargs.items()]
        print "%s(%s)" % (func.__name__, ", ".join(params))
        try:
            ret = func(*args, **kwargs)
        except Exception, e:
            ret = e
            raise
        finally:
            print "%s() -> %s" % (func.__name__, quote(ret))
        return ret

    return wrapper

def debugexc(func):
    """Decorator function. Prints the arguments and the return value of a function."""

    def quote(x):
        if type(x) == str: return '"%s"' % x[:100]
        else: return repr(x)[:300]

    def wrapper(*args, **kwargs):
        #traceback.print_tb()
        params = [quote(x) for x in args] + ["%s=%s" % (k, quote(v)) for k, v in kwargs.items()]
        try:
            ret = func(*args, **kwargs)
        except Exception, e:
            traceback.print_exc()
            ret = e
            raise
        return ret

    return wrapper

def _contr_str(container):
    if not container:
        return ''
    ret = '{'
    for k in container.keys():
        d,t = container[k]
        ret += (k + ': (')
        if len(d) > 128:
            ret += ('...data...len:' + str(len(d)) + '...')
        else:
            ret += str(d)
        ret += (', ' + t + ') ')
    ret += '}'
    return ret

def print_msg(op, subject, subj_id, arg_list, container):
    #print op, subject, subj_id, arg_list, container
    #print op, subject, subj_id, arg_list, _contr_str(container)
    pass

class DummyLock(object):
    def __init__(self, name=''):
        self.name = name
    def __enter__(self):
        print 'dummylock %s ENTER' % self.name
    def __exit__(self, *args):
        print 'dummylock %s EXIT' % self.name
