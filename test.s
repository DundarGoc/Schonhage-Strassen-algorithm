	.file	"test.c"
	.text
	.section	.rodata.str1.1,"aMS",@progbits,1
.LC1:
	.string	"%i %i %i %i %i %i %i %i\n"
	.section	.text.startup,"ax",@progbits
	.p2align 4
	.globl	main
	.type	main, @function
main:
.LFB5290:
	.cfi_startproc
	leaq	8(%rsp), %r10
	.cfi_def_cfa 10, 0
	andq	$-32, %rsp
	pushq	-8(%r10)
	xorl	%r9d, %r9d
	xorl	%r8d, %r8d
	pushq	%rbp
	xorl	%ecx, %ecx
	xorl	%edx, %edx
	.cfi_escape 0x10,0x6,0x2,0x76,0
	movq	%rsp, %rbp
	pushq	%r10
	.cfi_escape 0xf,0x3,0x76,0x78,0x6
	xorl	%esi, %esi
	movl	$.LC1, %edi
	subq	$80, %rsp
	pushq	$0
	xorl	%eax, %eax
	pushq	$0
	pushq	$0
	call	printf
	movq	-8(%rbp), %r10
	.cfi_def_cfa 10, 0
	addq	$32, %rsp
	xorl	%eax, %eax
	leave
	leaq	-8(%r10), %rsp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5290:
	.size	main, .-main
	.ident	"GCC: (GNU) 9.1.1 20190503 (Red Hat 9.1.1-1)"
	.section	.note.GNU-stack,"",@progbits
