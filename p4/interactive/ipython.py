from IPython import get_ipython

def custom_exc(shell, etype, evalue, tb, tb_offset=None):
    import traceback
    from p4.interactive.excepthook import invoke_editor 
    te = traceback.extract_tb(tb)
    shell.showtraceback((etype, evalue, tb), tb_offset=tb_offset)
    invoke_editor(te)

# Tell IPython to use it for any/all Exceptions:
get_ipython().set_custom_exc((Exception,), custom_exc)

