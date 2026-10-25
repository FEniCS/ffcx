	.file	"<stdin>"
	.text
	.p2align 4
	.globl	f
	.type	f, @function
f:
.LFB0:
	.cfi_startproc
	endbr64
	vfmadd132sd	%xmm1, %xmm2, %xmm0
	ret
	.cfi_endproc
.LFE0:
	.size	f, .-f
	.ident	"GCC: (Ubuntu 14.3.0-14ubuntu1) 14.3.0"
	.section	.note.GNU-stack,"",@progbits
	.section	.note.gnu.property,"a"
	.align 8
	.long	1f - 0f
	.long	4f - 1f
	.long	5
0:
	.string	"GNU"
1:
	.align 8
	.long	0xc0000002
	.long	3f - 2f
2:
	.long	0x3
3:
	.align 8
4:
