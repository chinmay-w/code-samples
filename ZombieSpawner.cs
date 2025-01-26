/*
    When working on horror-rhythm video game Zombeats (https://umbresp.itch.io/zombeats), I wanted to add the seeminly-minor 
    feature of a 4-bar count-off that would play before each of the songs in the game. This turned out to be incredibly annoying 
    for multiple reasons. Firstly, at the time our system for mapping out a level over time was tempo-blind (it was expressed 
    purely as a sequence of Zombies that spawn with a given time-offset from the start of the song, rather than abstracted 
    through beats and measures). 
    
    Secondly, Unity's handling of audio playback is actually quite convoluted, with many of the most common methods / properties 
    being unsuitable for anything requiring precision... like a count off, which needs to be on-beat. So in the process of this 
    seemingly minor addition I ended up heavily refactoring and implementing several global systems. But because they were 
    useful throughout the rest of development, I don't regret it--being able to specify audio timings more precisely helped a lot 
    with timing sound effects, and more precise tempo handling allowed us to implement features later (like being able to play 
    the same level at three different speeds) more straightforwardly.
*/

// (I've deleted all the code in the class unrelated to this particular feature)

public class ZombieSpawner : MonoBehaviour
{
    public const float GLOBAL_AUDIO_DELAY = 0.5f;
    // Unity has a really hard time playing audio with any kind of precision unless it has time beforehand to load the audio file
    // With 0.5 seconds of silence in which the engine can load everything, we were able to avoid a lot of weirdness on song start

    private const float COUNT_OFF_REPS = 4; // assumes all songs are in 4/4 time (initially)

    public bool songStarted;
    public float startTempo; // initial tempo of song
    private float currentTempo;
    public float Tempo
    {
        get => currentTempo;
    }

    private int firstNoteOfCurrentTempo; // index for beatmap.notes[], set internally

    // Our beatmap class tracks notes based on time from start of song
    // The getter finds the offset between song start and beatmap start (in seconds)
    public float FirstNoteOffset
    {
        // SongSpeedUtil keeps track of the 'speed modifier' that the player has selected for the current level
        get => (float)beatmap.notes[firstNoteOfCurrentTempo].time / SongSpeedUtil.GetSpeed(SongSpeedUtil.CurrentSongSpeed);
    }

    private AudioSource[] drumstickSfx; // the audio sources for the four drumstick sfx involved in the count-off
    private AudioSource _musicPlayer; // Audio source for the level's song
    private Beatmap beatmap; // Our abstraction of the level's beatmap once loaded from json

    void Awake()
    {
        songStarted = false;
        drumstickSfx = transform.GetChild(0).GetComponents<AudioSource>();
    }

    IEnumerator trackSongStart(float delay) // coroutine
    {
        yield return new WaitForSeconds(delay); 
        songStarted = true;
    }

    // We cannot use Unity's Start OR Awake functions to implement musical logic because of the loading stutter
    // Instead, we created a custom third pass once loading is complete (before that we throw up a loading screen)
    // In 'LoadingBlocker.cs': GameObject.Find("Level").BroadcastMessage("InitAfterLoad"); 
    // which runs InitAfterLoad() on all objects synchronously (if it exists)
    void InitAfterLoad()
    {
        beatmap = Beatmap.fromJson(_beatmapJson.text);
        beatmapLength = beatmap.notes.Count;

        _musicPlayer.volume = PlayerPrefs.GetFloat("musicvol", 0.5f);

        startTempo = PlayerPrefs.GetInt("startTempo"); 
        // PlayerPrefs are used to communicate between scenes
        // In this case, the start tempo is directly shown to the player on the level select
        // and then retrieved on level-load (where we are currently)

        startTempo *= SongSpeedUtil.GetSpeed(SongSpeedUtil.CurrentSongSpeed);
        // An example of a feature we added later using the systems I made for this minor feature
        // We modify the tempo based on the player's difficulty settings and can spawn Zombies accordingly later
        currentTempo = startTempo;
        firstNoteOfCurrentTempo = 0;

        int beatLength = 60 / Tempo;
        // length of one beat in seconds, e.g. for a song with a BPM of 120 each beat would take 0.5 seconds
        songStarted = false;

        if (FirstNoteOffset < COUNT_OFF_REPS * beatLength) // handles the case where the beatmap starts < 4 beats into the song
        {
            // AudioSource.PlayScheduled() plays songs at a far more precise time than either Play() or PlayDelayed()
            // AudioSettings.dspTime is a far more accurate timer than Time.time
            for (int i = 0; i < COUNT_OFF_REPS; i++) {
                drumstickSfx[i].PlayScheduled(AudioSettings.dspTime + i * beatLength + GLOBAL_AUDIO_DELAY);
            }

            _musicPlayer.PlayScheduled(AudioSettings.dspTime + (COUNT_OFF_REPS * beatLength) - FirstNoteOffset + GLOBAL_AUDIO_DELAY);
            StartCoroutine(trackSongStart((4 * beatLength) - FirstNoteOffset + GLOBAL_AUDIO_DELAY));
            // starts the song such that the next bar downbeat of the count-off is the first note in the beatmap
        }
        else // if the beatmap starts later than 4 beats in, then we delay the count-off until it's about to start
        {
            _musicPlayer.PlayScheduled(AudioSettings.dspTime + GLOBAL_AUDIO_DELAY);
            StartCoroutine(trackSongStart(GLOBAL_AUDIO_DELAY));
            // start song asap

            for (int i = 1; i <= COUNT_OFF_REPS; i++)
            {
                drumstickSfx[i].PlayScheduled(AudioSettings.dspTime + FirstNoteOffset - (i * beatLength) + GLOBAL_AUDIO_DELAY);
            }
            // chedule count-off for the four bars before the first note plays
        }
    }
}