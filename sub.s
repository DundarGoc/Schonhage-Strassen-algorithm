	.file	"sub.c"
	.text
	.p2align 4
	.globl	test
	.type	test, @function
test:
.LFB5279:
	.cfi_startproc
	vmovdqu	(%rsi), %ymm0
	vpsubd	(%rdx), %ymm0, %ymm0
	vmovdqu	%ymm0, (%rdi)
	vzeroupper
	ret
	.cfi_endproc
.LFE5279:
	.size	test, .-test
	.ident	"GCC: (GNU) 9.1.0"
	.section	.note.GNU-stack,"",@progbits
