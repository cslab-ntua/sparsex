#include <linux/module.h>
#include <linux/kernel.h>
#include <asm/processor.h>

/*
 * simple module to enable PCE flag in CR4 register
 * this flag will allow to use the rdpmc instruction
 * to read perfomance counters
 *
 * --kkourt@cslab.ece.ntua.gr
 */

static void cr4_pce_enable(void *arg)
{
	set_in_cr4(X86_CR4_PCE);
}

static void cr4_pce_disable(void *arg)
{
	clear_in_cr4(X86_CR4_PCE);
}

int cr4_pce_enable_init(void)
{

	printk(KERN_ALERT "cr4_pce_enable: enabling pce in cr4\n");

	#ifdef CONFIG_X86
	on_each_cpu(cr4_pce_enable, NULL, 0, 0);
	#else
	print(KERN_ALRT "FAILED: not x86\n");
	#endif
		
	return 0;
}

void cr4_pce_enable_exit(void)
{
	printk(KERN_ALERT "cr4_pce_enable: disabling pce in cr4\n");

	#ifdef CONFIG_X86
	on_each_cpu(cr4_pce_disable, NULL, 0, 0);
	#else
	print(KERN_ALRT "FAILED: not x86\n");
	#endif
}

module_init(cr4_pce_enable_init);
module_exit(cr4_pce_enable_exit);
MODULE_LICENSE("GPL");
