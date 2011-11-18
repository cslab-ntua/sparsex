from datetime import *
from re import *
from string import *
from sys import *

global program_name

available_operations = {
    'generate_bcsr' : 'template, r, c, outfile',
    }

def print_available_operations(fout):
    fout.write('Available operations:\n')
    for name,args in available_operations.iteritems():
        fout.write('\t%s(%s)\n' % (name, args))

def print_usage(fout):
    global program_name
    fout.write('Usage: %s <OP> [option=value] ...\n' % program_name)
    print_available_operations(fout)

def generate_bcsr(template, r, c, outfile):
    global program_name
    r = int(r)
    c = int(c)
    if r <= 0 or c <=0:
        raise UserWarning('invalid block dimensions')
    
    f_template = open(template, 'r')
    template_text = ''
    for line in f_template:
        template_text = template_text + line

    src_template = Template(template_text)

    # generate the loop body
    y_init = '\t\t'
    loop_body = '\t\t\t'
    y_assign = '\t\t'
    for i in xrange(0,r):
        y_init = y_init + 'register ELEM_TYPE yr%d = 0;' % i
        if i < r-1:
            y_init = y_init + '\n\t\t'
        loop_body = loop_body + 'yr%d += ' % i
        for j in xrange(0,c):
            vi = i*c + j
            loop_body = loop_body + 'bvalues[j+%d]*x[x_start+%d]' % (vi, j)
            if j == c-1:
                loop_body = loop_body + ';'
                if i < r-1:
                    loop_body = loop_body + '\n\t\t\t'
            else:
                loop_body = loop_body + ' + ';
        y_assign = y_assign + 'y[i+%d] = yr%d;' % (i, i)
        if i < r-1:
           y_assign = y_assign + '\n\t\t'

    src_out = src_template.substitute(current_year = date.today().year,
                                      generator_prog = program_name,
                                      src_filename = outfile,
                                      r = r, c = c,
                                      y_init_hook = y_init,
                                      loop_body_hook = loop_body,
                                      y_assign_hook = y_assign);
    # do some aesthetic postprocessing
    src_out = sub('\+0', '', src_out)
    src_out = sub('\s+\/\s+1', '', src_out)
    src_out = sub('\t', '    ', src_out)

    # write output
    if outfile == 'stdout':
        fout = stdout
    elif outfile == 'stderr':
        fout = stderr
    else:
        fout = open(outfile, 'w')
    fout.write(src_out)

if __name__ == '__main__':
    global program_name
    argc = len(argv)
    program_name = argv[0]
    if argc < 2:
        stderr.write('%s: too few arguments\n' % program_name)
        print_usage(stderr)
        exit(1)

    op = argv[1]
    op_args = argv[2:]
    op_kwargs = dict()
    for arg in op_args:
        key, val = arg.split('=', 2)
        op_kwargs[key] = val

    if not op in available_operations.keys():
        stderr.write('unknown operation: %s\n' % op)
        print_available_operations(stderr)
        exit(1)

    locals()[op](**op_kwargs)
