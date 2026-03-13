# Задание 1: Парсинг FASTQ файла
def parse_fastq(filename):
    records = []
    with open(filename, 'r') as f:
        lines = f.readlines()
    for i in range(0, len(lines), 4):
        if i + 3 >= len(lines):
            break
        identifier = lines[i].strip()[1:].split()[0]
        sequence = lines[i + 1].strip()
        quality = lines[i + 3].strip()

        records.append({
            'id': identifier,
            'sequence': sequence,
            'quality': quality,
            'length': len(sequence)
        })
    return records

# Задание 1: Статистика по прочтениям
def get_statistics(records):
    lengths = [r['length'] for r in records]
    return {
        'total_reads': len(records),
        'min_length': min(lengths),
        'max_length': max(lengths),
        'avg_length': round(sum(lengths) / len(lengths))
    }

# Задание 2: Расчет GC-состава
def calculate_gc_content(records):
    total_gc = 0
    total_bases = 0

    for record in records:
        sequence = record['sequence'].upper()
        total_gc += sequence.count('G') + sequence.count('C')
        total_bases += len(sequence)
    return round((total_gc / total_bases) * 100, 2) if total_bases > 0 else 0.0
# Задание 3: Среднее качество Phred для позиции
def calculate_avg_phred_position(records, position=10):
    phred_scores = []

    for record in records:
        quality = record['quality']
        if len(quality) >= position:
            ascii_char = quality[position - 1]
            phred_score = ord(ascii_char) - 33
            phred_scores.append(phred_score)

    return round(sum(phred_scores) / len(phred_scores)) if phred_scores else 0

# Задание 4: Применение тримминга по качеству
def apply_quality_trimming(records, window_size=5, min_quality=30):
    trimmed_records = []
    trimmed_count = 0
    for record in records:
        quals = [ord(char) - 33 for char in record['quality']]
        keep_length = None
        if len(quals) >= window_size:
            threshold = min_quality * window_size
            window_sum = sum(quals[:window_size])

            if window_sum >= threshold:
                keep_length = len(quals)
                for i in range(len(quals) - window_size):
                    window_sum += quals[i + window_size] - quals[i]
                    if window_sum < threshold:
                        keep_length = i + window_size
                        break
                while keep_length > 1 and quals[keep_length - 1] < min_quality:
                    keep_length -= 1
                if keep_length <= 0:
                    keep_length = None

        if keep_length is not None and keep_length > 0:
            trimmed_records.append({
                'id': record['id'],
                'sequence': record['sequence'][:keep_length],
                'quality': record['quality'][:keep_length],
                'length': keep_length,
                'original_length': record['length']
            })
            if keep_length < record['length']:
                trimmed_count += 1
    return trimmed_records, trimmed_count

# Задание 4: Фильтрация по длине
def filter_by_length(records, min_length=60):
    return [r for r in records if r['length'] >= min_length]

filename = "reads.fastq.txt"
# Вывод задания 1: Парсинг и статистика
records = parse_fastq(filename)
stats = get_statistics(records)

print(f"Общее число прочтений в файле равно {stats['total_reads']}.")
print(f"Минимальная длина прочтения равна {stats['min_length']}.")
print(f"Средняя длина прочтения равна {stats['avg_length']}.")
print(f"Максимальная длина прочтения равна {stats['max_length']}.")

# Вывод задания 2: GC-состав
gc_content = calculate_gc_content(records)
print(f"GC-состав: {gc_content}")

# Вывод задания 3: Среднее качество Phred для позиции 10
avg_phred = calculate_avg_phred_position(records, 10)
print(f"Среднее качество Phred для позиции 10: {avg_phred}")

# Вывод задания 4: Тримминг по качеству (окно 5, качество 30)
trimmed_records, trimmed_count = apply_quality_trimming(records, window_size=5, min_quality=30)
print(f"\nПосле тримминга по качеству:")
print(f"Прочтений подверглось триммингу: {trimmed_count}")
trimmed_stats = get_statistics(trimmed_records)
print(f"Минимальная длина: {trimmed_stats['min_length']}")
print(f"Средняя длина: {trimmed_stats['avg_length']}")
print(f"Максимальная длина: {trimmed_stats['max_length']}")
#Вывод задания 4: Фильтрация по длине >= 60
filtered_records = filter_by_length(trimmed_records, min_length=60)
print(f"\nПосле фильтрации по длине >= 60:")
print(f"Осталось прочтений: {len(filtered_records)}")
