	.file	"sse.c"
	.text
	.globl	sse_add
	.type	sse_add, @function
sse_add:
.LFB653:
	.cfi_startproc
	endbr64
	testq	%rcx, %rcx
	je	.L1
	movl	$0, %eax
.L3:
	movaps	(%rdi,%rax,4), %xmm0
	addps	(%rsi,%rax,4), %xmm0
	movaps	%xmm0, (%rdx,%rax,4)
	addq	$4, %rax
	cmpq	%rax, %rcx
	ja	.L3
.L1:
	ret
	.cfi_endproc
.LFE653:
	.size	sse_add, .-sse_add
	.globl	sse_add_d
	.type	sse_add_d, @function
sse_add_d:
.LFB654:
	.cfi_startproc
	endbr64
	testq	%rcx, %rcx
	je	.L5
	movl	$0, %eax
.L7:
	movapd	(%rdi,%rax,8), %xmm0
	addpd	(%rsi,%rax,8), %xmm0
	movaps	%xmm0, (%rdx,%rax,8)
	addq	$2, %rax
	cmpq	%rax, %rcx
	ja	.L7
.L5:
	ret
	.cfi_endproc
.LFE654:
	.size	sse_add_d, .-sse_add_d
	.ident	"GCC: (Ubuntu 11.3.0-1ubuntu1~22.04) 11.3.0"
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
