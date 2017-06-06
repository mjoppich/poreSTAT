

def calcN50(lengths):
    slengths = sorted(lengths, reverse=True)
    totalLength = sum(slengths)

    currentLength = 0
    neededLength = totalLength / 2.0

    print(neededLength)

    n50Value = 0
    L50 = 0

    for i in range(0, len(slengths)):

        x = slengths[i]

        currentLength += x
        L50 += 1

        if currentLength >= neededLength:
            n50Value = x

            print("N50 value is " + str(x) + " for length " + str(totalLength) + " (in half: " + str(
                neededLength) + " ) and L50: " + str(i) + " " + str(len(slengths)))
            break

    return (n50Value, L50)