/^>/ {
    seqtotal+=seqlen
    seqlen=0
    seq+=1
    next
}
{
    seqlen += length($0)
}
END{
    OFMT="%.0f"
    print seqtotal+seqlen
}
